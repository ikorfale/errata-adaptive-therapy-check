"""CM-CANCER-Q01 bench v3: aria's independent implementation of errata's v3 spec (Abund a76c9913, 2026-10-01;
spec owner errata, PLAN #70). Written from the spec text; errata's code was read, not run.

Model (Zhang 2017 / Strobl 2021 form, PSA units), daily steps in log space:
  dS = rS S g (1 - dD D) - dT S,   dR = rR R g - dT R,   g = 1 - (S+R)/K,   burden N = S+R.
Parameter sets: every fit inside a patient's 95% profile band (errata's profile_all.json, tag v1.0 / ad66924,
field p = log [rS, rR, dT, dD, K, S0, R0]); 67 patients of the Bruchovsky 2006 IAS trial.
Rules decide at 7-day visits from burden only (x = N/N0); the dose is held between visits.
Progression: N > 1.2 N0. Horizon H = 1825 d.
Score on one set: an entry that progresses before H FAILS. Else score = frontier dose / entry dose, where the frontier
is the least cumulative dose of any containment rule (weekly modulation, targets 0.10..1.18 step 0.02, gain 2;
plus on/off: off at x <= lo, on at x >= hi) that survives H at the same or lower mean burden.
VACUOUS set: the untreated tumour never crosses the line within H (any rule survives); dropped and counted.
An entry wins a patient only if score > 1.01 on EVERY non-vacuous set of that patient. Report patients, not means.
Usage: python3 results/cm_cancer_q01_v3.py   (numpy; writes results/CM-CANCER-Q01-v3.json)"""
import csv, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "cancer_q01")
H, VISIT, LINE, WIN = 1825, 7, 1.2, 1.01
TARGETS = [round(0.10 + 0.02 * i, 2) for i in range(55)]                       # 0.10 .. 1.18
ONOFF = [(round(lo, 1), round(hi, 1)) for lo in np.arange(0.1, 1.01, 0.1) for hi in np.arange(0.2, 1.16, 0.1) if hi > lo + 0.05]

def trial():
    """Patient on/off series from the public file (columns: id, date, CPA, LEU, PSA, T, cycle, on, day, day2)."""
    out = {}
    for f in sorted(glob.glob(os.path.join(DATA, "bruchovsky", "dataTanaka", "Bruchovsky_et_al", "patient*.txt"))):
        rows = []
        for r in csv.reader(open(f)):
            try: rows.append((float(r[8].strip('"')), float(r[7].strip('"'))))
            except (ValueError, IndexError): pass
        rows.sort(); d0 = rows[0][0]
        out[os.path.basename(f)[7:10]] = [(d - d0, on) for d, on in rows]
    return out

def first_cycle(rows):
    """The patient's own first on phase and first off phase (days); off phase 365 d if treatment never resumed."""
    off = next(t for t, on in rows if on == 0)
    back = next((t for t, on in rows if on == 1 and t > off), None)
    return off, (back - off) if back is not None else 365.0

# entries: name -> rule(t, x, state) -> dose, x = N/N0 at a visit
def mtd(t, x, s): return 1.0
def none(t, x, s): return 0.0
def adaptive50(t, x, s):
    if s.get("on", True) and x <= 0.5: s["on"] = False
    elif not s.get("on", True) and x >= 1.0: s["on"] = True
    return float(s.get("on", True))
def modulate(target):
    def rule(t, x, s):
        s["D"] = min(1.0, max(0.0, s.get("D", 1.0) + 2.0 * (x - target))); return s["D"]
    return rule
def onoff(lo, hi):
    def rule(t, x, s):
        if s.get("on", True) and x <= lo: s["on"] = False
        elif not s.get("on", True) and x >= hi: s["on"] = True
        return float(s.get("on", True))
    return rule
def setpoint_088(t, x, s):
    """agentcue CM-CANCER-103 in bench-v2 units (0.88 absolute at N0 0.75 = 1.1733 N0; gain 20 per 1.0 = 15 per N0)."""
    return min(1.0, max(0.0, 15.0 * (x - 0.88 / 0.75)))
def replay_cycle1(c1):
    on_d, off_d = c1
    def rule(t, x, s): return float(t % (on_d + off_d) < on_d)
    return rule

def simulate(p, rule):
    """Returns (ttp, cumulative dose, mean burden in N0) up to progression or H."""
    rS, rR, dT, dD, K, S0, R0 = np.exp(p)
    lS, lR, N0, D, state, dose, area = np.log(S0), np.log(R0), S0 + R0, 1.0, {}, 0.0, 0.0
    for t in range(H):
        S, R = np.exp(lS), np.exp(lR); N = S + R
        if N > LINE * N0: return t, dose, area / max(t, 1)
        if t % VISIT == 0: D = rule(t, N / N0, state)
        dose += D; area += N / N0
        g = 1.0 - N / K
        lS = max(lS + rS * g * (1.0 - dD * D) - dT, -40.0); lR = max(lR + rR * g - dT, -40.0)
    return H, dose, area / H

def score_set(p, entries):
    fam = [simulate(p, modulate(x)) for x in TARGETS] + [simulate(p, onoff(lo, hi)) for lo, hi in ONOFF]
    fam = [(mb, d) for ttp, d, mb in fam if ttp >= H]
    vacuous = simulate(p, none)[0] >= H
    out = {}
    for name, rule in entries.items():
        ttp, d, mb = simulate(p, rule)
        if ttp < H: out[name] = ("FAIL", ttp); continue
        ok = [fd for fmb, fd in fam if fmb <= mb + 1e-9]
        out[name] = ("OK", min(ok) / max(d, 1e-9), round(mb, 3), round(d, 1)) if ok else ("BELOW", None, round(mb, 3), round(d, 1))
    return vacuous, out

if __name__ == "__main__":
    prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
    pts = trial()
    names = ["mtd", "adaptive50", "replay_cycle1", "mod_1.00", "mod_0.50", "setpoint_088"]
    res, tally = {}, {n: dict(win=0, fail_some=0, below_some=0, scored=0) for n in names}
    n_sets = n_vac = 0; usable_all = usable_any = 0
    for pid in sorted(prof):
        band = [q["p"] for q in prof[pid]["pts"] if q["ok"]]
        entries = {"mtd": mtd, "adaptive50": adaptive50, "replay_cycle1": replay_cycle1(first_cycle(pts[pid])),
                   "mod_1.00": modulate(1.0), "mod_0.50": modulate(0.5), "setpoint_088": setpoint_088}
        sets = [score_set(np.array(p), entries) for p in band]
        n_sets += len(sets); vac = [v for v, _ in sets]; n_vac += sum(vac)
        usable_all += not any(vac); usable_any += not all(vac)
        live = [o for v, o in sets if not v]
        res[pid] = dict(n_band=len(sets), n_vacuous=sum(vac), sets=[{k: list(x) for k, x in o.items()} for _, o in sets])
        if not live: continue
        for n in names:
            sc = [o[n] for o in live]; tally[n]["scored"] += 1
            tally[n]["fail_some"] += any(s[0] == "FAIL" for s in sc)
            tally[n]["below_some"] += any(s[0] == "BELOW" for s in sc)
            tally[n]["win"] += all(s[0] == "OK" and s[1] > WIN for s in sc)
    summary = dict(patients=len(prof), band_sets=n_sets, vacuous_sets=n_vac,
                   patients_no_vacuous_set=usable_all, patients_some_live_set=usable_any, entries=tally)
    json.dump(dict(spec="errata v3 (Abund a76c9913, repo ad66924)", summary=summary, patients=res),
              open(os.path.join(HERE, "CM-CANCER-Q01-v3.json"), "w"), indent=0)
    print(json.dumps(summary, indent=1))
