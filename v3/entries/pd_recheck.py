"""Re-check the three robust wins in pd_entry.json against a finer frontier, the way 106 was re-checked:
on-off hysteresis on a 0.02 grid and integral/proportional with more gains, each kept only if it contains the
tumour for H at mean burden <= the entry's. Score = best family dose / entry dose (> 1.01 = entry still wins)."""
import json, os, sys, numpy as np
from cm_cancer_q01_v3 import DATA, simulate, onoff, WIN
from cm_cancer_q01_v31 import integral, proportional
from pd_entry import ENTRIES
from cm_cancer_q01_v3 import H
prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
none = lambda t, x, s: 0.0
CASES = [("064", "trend_onoff_0.7_1.10"), ("088", "pid_0.5_d10"), ("064", "trend_onoff_0.9_1.15")]
grid = np.round(np.arange(0.30, 1.19, 0.02), 2)
fine_onoff = [(lo, hi) for lo in grid for hi in grid if hi > lo + 0.01]
for pid, name in CASES:
    f, a = ENTRIES[name]
    for k, p in enumerate(q["p"] for q in prof[pid]["pts"] if q["ok"]):
        p = np.array(p)
        if simulate(p, none)[0] >= H: continue                         # vacuous fit, not scored
        ttp, d, mb = simulate(p, f(*a))
        if ttp < H: print(pid, name, "fit", k, "FAIL"); continue
        best = None; which = None
        fam = [("onoff", lh, onoff(*lh)) for lh in fine_onoff] + \
              [("int", (t, g), integral(t, g)) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in np.round(np.arange(mb - 0.2, mb + 0.01, 0.01), 3)] + \
              [("prop", (t, g), proportional(t, g)) for g in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in np.round(np.arange(mb - 0.3, mb + 0.01, 0.01), 3)]
        for kind, par, r in fam:
            t2, d2, mb2 = simulate(p, r)
            if t2 >= H and mb2 <= mb + 1e-9 and (best is None or d2 < best): best, which = d2, (kind, par)
        print(pid, name, "fit", k, "entry dose %.1f mb %.4f | finer frontier dose %s by %s | score %s" %
              (d, mb, None if best is None else round(best, 1), which, None if best is None else round(best / d, 4)), flush=True)
