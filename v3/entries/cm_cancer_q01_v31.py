"""CM-CANCER-Q01 bench v3.1 (errata, 2026-10-02): v3 plus the two rulings on aria's holes (Abund 08e598b4 -> aba281df,
a63c3170). Model, parameter sets, entries, line and horizon are imported unchanged from cm_cancer_q01_v3.py.

Changes against v3:
1. Refined frontier. The containment family adds integral modulation D += g (x - target) with gains 2/5/15 and
   proportional setpoints D = clamp(g (x - target)) with gains 5/15/30, targets every 0.002 from 1.000 to 1.198.
2. Local refinement. Any entry that still scores > 1.01 on a set is re-scored against a grid around its own mean
   burden mb: targets mb-0.01 .. mb+0.01 every 0.0002, integral gains 2/5/15, proportional gains 5/10/15/20/30.
   A score that falls to <= 1.01 there was a grid artefact. The reported score is the lowest of the three.
   v3.2 (2026-10-03, after trend-reacting entries showed false wins): local refinement also re-tests on-off
   hysteresis with lo/hi on a 0.02 grid (lo 0.30-1.18, hi > lo), and widens the gains to integral 1/2/5/10/15/30
   and proportional 3/5/10/15/20/30/60. The three v3.1 wins by trend entries all fall below 1.01 under it.
3. Live sets. A set is live if the untreated tumour crosses the line within H. Vacuous sets rank nothing, so their
   score is null (dividing by a near-zero dose gave values like 5e9). A patient is eligible with >= 2 live sets.
   An entry wins an eligible patient only if it scores > 1.01 on EVERY live set.
Usage: python3 results/cm_cancer_q01_v31.py   (numpy; writes results/CM-CANCER-Q01-v31.json; about 3 min on one CPU)"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from cm_cancer_q01_v3 import (DATA, H, WIN, TARGETS, ONOFF, trial, first_cycle, simulate, mtd, none, adaptive50,
                              modulate, onoff, setpoint_088, replay_cycle1)

MIN_LIVE = 2
_G = [round(0.30 + 0.02 * i, 2) for i in range(45)]                             # 0.30 .. 1.18
ONOFF_FINE = [(lo, hi) for lo in _G for hi in _G if hi > lo + 0.01]             # v3.2: 990 on-off pairs
FINE = [round(1.000 + 0.002 * i, 3) for i in range(100)]                       # 1.000 .. 1.198

def integral(target, g):
    def rule(t, x, s):
        s["D"] = min(1.0, max(0.0, s.get("D", 1.0) + g * (x - target))); return s["D"]
    return rule
def proportional(target, g):
    def rule(t, x, s): return min(1.0, max(0.0, g * (x - target)))
    return rule

def best(runs, mb):
    """Least dose among simulated family runs (ttp, dose, mean burden) that survive H at mean burden <= mb."""
    ok = [d for ttp, d, fmb in runs if ttp >= H and fmb <= mb + 1e-9]
    return min(ok) if ok else None

def score_set(p, entries, refined):
    if simulate(p, none)[0] >= H:
        return False, {n: ["VACUOUS", None] for n in entries}
    fam = [simulate(p, r) for r in refined]                                      # v3 family first, then the fine one
    coarse = fam[:len(TARGETS) + len(ONOFF)]
    out = {}
    for name, rule in entries.items():
        ttp, d, mb = simulate(p, rule)
        if ttp < H: out[name] = ["FAIL", ttp]; continue
        sc = {}
        for tag, runs in (("v3", coarse), ("refined", fam)):
            b = best(runs, mb); sc[tag] = round(b / max(d, 1e-9), 4) if b is not None else None
        if (sc["refined"] or 0) > WIN:
            ts = np.round(np.arange(mb - 0.01, mb + 0.0101, 0.0002), 5)
            local = [integral(t, g) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in ts] + \
                    [proportional(t, g) for g in (3.0, 5.0, 10.0, 15.0, 20.0, 30.0, 60.0) for t in ts] + \
                    [onoff(lo, hi) for lo, hi in ONOFF_FINE]
            b = best([simulate(p, r) for r in local], mb); sc["local"] = round(b / max(d, 1e-9), 4) if b is not None else None
        vals = [v for v in sc.values() if v is not None]
        out[name] = ["OK", min(vals) if vals else None, round(mb, 4), round(d, 1), sc]
    return True, out

if __name__ == "__main__":
    prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
    pts = trial()
    refined = [modulate(x) for x in TARGETS] + [onoff(lo, hi) for lo, hi in ONOFF] + \
              [integral(t, g) for g in (2.0, 5.0, 15.0) for t in FINE] + \
              [proportional(t, g) for g in (5.0, 15.0, 30.0) for t in FINE]
    names = ["mtd", "adaptive50", "replay_cycle1", "mod_1.00", "mod_0.50", "setpoint_088"]
    tally = {n: dict(win=0, fail_some=0, below_some=0) for n in names}
    res, live_hist = {}, {}
    for pid in sorted(prof):
        band = [q["p"] for q in prof[pid]["pts"] if q["ok"]]
        entries = {"mtd": mtd, "adaptive50": adaptive50, "replay_cycle1": replay_cycle1(first_cycle(pts[pid])),
                   "mod_1.00": modulate(1.0), "mod_0.50": modulate(0.5), "setpoint_088": setpoint_088}
        sets = [score_set(np.array(p), entries, refined) for p in band]
        live = [o for v, o in sets if v]
        live_hist[len(live)] = live_hist.get(len(live), 0) + 1
        res[pid] = dict(n_band=len(sets), n_live=len(live), eligible=len(live) >= MIN_LIVE, sets=[o for _, o in sets])
        print(pid, len(sets), "live", len(live), flush=True)
        if len(live) < MIN_LIVE: continue
        for n in names:
            sc = [o[n] for o in live]
            tally[n]["fail_some"] += any(s[0] == "FAIL" for s in sc)
            tally[n]["below_some"] += any(s[0] == "OK" and s[1] is None for s in sc)
            tally[n]["win"] += all(s[0] == "OK" and s[1] is not None and s[1] > WIN for s in sc)
    summary = dict(patients=len(prof), patients_by_live_sets=dict(sorted(live_hist.items())),
                   eligible=sum(r["eligible"] for r in res.values()), min_live_sets=MIN_LIVE, entries=tally)
    json.dump(dict(spec="errata v3.2 (Abund aba281df, a63c3170; on-off local refinement)", summary=summary, patients=res),
              open(os.path.join(HERE, "CM-CANCER-Q01-v31.json"), "w"), indent=0)
    print(json.dumps(summary, indent=1))
