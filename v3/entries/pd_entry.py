"""Entry test for board post 71589: does a rule that reacts to the TREND of burden (the one shape the v3.1 frontier
family lacks) beat the frontier on any eligible patient? Scored with the unchanged v3.1 harness (score_set)."""
import json, os, sys, numpy as np
from cm_cancer_q01_v3 import DATA, TARGETS, ONOFF, WIN, modulate, onoff
from cm_cancer_q01_v31 import MIN_LIVE, FINE, integral, proportional, score_set

def pid_rule(target, gi, gd):
    """D += gi*(x - target) + gd*(x - x_prev): integral on level plus a kick on the change since the last visit."""
    def rule(t, x, s):
        xp = s.get("x", x); s["x"] = x
        s["D"] = min(1.0, max(0.0, s.get("D", 1.0) + gi * (x - target) + gd * (x - xp))); return s["D"]
    return rule
def trend_onoff(lo, hi):
    """On-off with hysteresis, but switch on early when burden is rising and above lo (trend-aware)."""
    def rule(t, x, s):
        xp = s.get("x", x); s["x"] = x; on = s.get("on", True)
        if on and x <= lo: on = False
        elif not on and (x >= hi or (x > lo and x - xp > 0.03)): on = True
        s["on"] = on; return float(on)
    return rule

ENTRIES = {}
for tg in (0.5, 0.8, 1.0, 1.1):
    for gd in (2.0, 5.0, 10.0):
        ENTRIES["pid_%.1f_d%g" % (tg, gd)] = (pid_rule, (tg, 2.0, gd))
for lo, hi in ((0.5, 1.0), (0.7, 1.1), (0.9, 1.15)):
    ENTRIES["trend_onoff_%.1f_%.2f" % (lo, hi)] = (trend_onoff, (lo, hi))

if __name__ == "__main__":
    elig = [k for k, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "CM-CANCER-Q01-v31.json")))["patients"].items() if v["eligible"]]
    prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
    refined = [modulate(x) for x in TARGETS] + [onoff(lo, hi) for lo, hi in ONOFF] + \
              [integral(t, g) for g in (2.0, 5.0, 15.0) for t in FINE] + \
              [proportional(t, g) for g in (5.0, 15.0, 30.0) for t in FINE]
    tally = {n: dict(win=0, fail_some=0, below_some=0, best_min=None) for n in ENTRIES}
    res = {}
    for pid in elig:
        band = [q["p"] for q in prof[pid]["pts"] if q["ok"]]
        sets = [score_set(np.array(p), {n: f(*a) for n, (f, a) in ENTRIES.items()}, refined) for p in band]
        live = [o for v, o in sets if v]; res[pid] = live
        for n in ENTRIES:
            sc = [o[n] for o in live]
            tally[n]["fail_some"] += any(s[0] == "FAIL" for s in sc)
            tally[n]["below_some"] += any(s[0] == "OK" and s[1] is None for s in sc)
            w = all(s[0] == "OK" and s[1] is not None and s[1] > WIN for s in sc); tally[n]["win"] += w
            if all(s[0] == "OK" and s[1] is not None for s in sc):
                m = min(s[1] for s in sc); b = tally[n]["best_min"]
                if b is None or m > b[0]: tally[n]["best_min"] = (m, pid)
        print(pid, "live", len(live), {n: [s[n][0] if s[n][0] != "OK" else s[n][1] for s in live] for n in list(ENTRIES)[:3]}, flush=True)
    json.dump(dict(tally=tally, patients=res), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "pd_entry.json"), "w"), indent=0)
    for n, t in tally.items(): print("%-22s win %d | fail on some fit %2d | below range %2d | best min-score %s" % (n, t["win"], t["fail_some"], t["below_some"], t["best_min"]))
