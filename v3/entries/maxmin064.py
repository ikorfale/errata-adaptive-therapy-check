"""zenith-claude 71607: give the trend on-off entry the same 0.02 refinement the frontier got. On patient 064,
search (lo, hi) on the 0.02 grid for the pair that maximises the MINIMUM score over the live fits, where
score = finer-frontier dose / entry dose at the entry's mean burden (same finer frontier as pd_recheck.py:
on-off 0.02 grid, integral and proportional setpoints on a 0.01 grid with the gains used there).
In-sample: the pair is tuned on the very fits it is scored on."""
import json, os, sys, numpy as np
from cm_cancer_q01_v3 import DATA, simulate, onoff, H, WIN
from cm_cancer_q01_v31 import integral, proportional
from pd_entry import trend_onoff
pid = sys.argv[1] if len(sys.argv) > 1 else "064"
prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
none = lambda t, x, s: 0.0
fits = [np.array(q["p"]) for q in prof[pid]["pts"] if q["ok"]]
fits = [p for p in fits if simulate(p, none)[0] < H]                    # live fits only
grid = np.round(np.arange(0.30, 1.19, 0.02), 2)
pairs = [(lo, hi) for lo in grid for hi in grid if hi > lo + 0.01]
sp = np.round(np.arange(0.30, 1.21, 0.01), 2)
fam = [onoff(*lh) for lh in pairs] + [integral(t, g) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in sp] + \
      [proportional(t, g) for g in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in sp]
F = []                                                                  # per fit: contained family runs (mb, dose)
for p in fits:
    runs = [simulate(p, r) for r in fam]
    F.append(np.array([(mb, d) for ttp, d, mb in runs if ttp >= H]))
def frontier(k, mb):
    ok = F[k][F[k][:, 0] <= mb + 1e-9]
    return ok[:, 1].min() if len(ok) else None
print(pid, "live fits", len(fits), "family rules", len(fam), "contained per fit", [len(f) for f in F], flush=True)
res = []
for lo, hi in pairs:
    sc = []
    for k, p in enumerate(fits):
        ttp, d, mb = simulate(p, trend_onoff(lo, hi))
        if ttp < H: sc.append(0.0); continue                            # tumour crosses: worst possible
        b = frontier(k, mb); sc.append(np.inf if b is None else b / d)
    res.append((min(sc), lo, hi, sc))
res.sort(key=lambda r: -r[0])
print("top 10 by min score over live fits:")
for m, lo, hi, sc in res[:10]: print(f"  lo {lo:.2f} hi {hi:.2f}  min {m:.4f}  per fit {[round(float(s), 4) for s in sc]}")
for k in range(len(fits)):
    r = max(res, key=lambda r: r[3][k])
    print(f"best pair for fit {k}: lo {r[1]:.2f} hi {r[2]:.2f}  scores on all fits {[round(float(s), 4) for s in r[3]]}")
print("pairs that contain the tumour on all live fits:", sum(1 for r in res if r[0] > 0), "of", len(res))
print("verdict:", "ROBUST WIN (in-sample)" if res[0][0] > WIN else "no robust win; max-min %.4f" % res[0][0])
