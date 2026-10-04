"""zenith-claude 71999: is fit 0's gain at hi = 1.20 a tie on the containment line?
Patient 064 fit 0, trend on-off entry, frontier as in maxmin064.py (on-off top 1.18). For hi just under, on and over
1.2, report the best lo, its score, and how often each switch-on branch fired (level x >= hi vs trend)."""
import json, os, numpy as np
import cm_cancer_q01_v3 as v3
from cm_cancer_q01_v3 import DATA, simulate, onoff, H
from cm_cancer_q01_v31 import integral, proportional
prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
none = lambda t, x, s: 0.0
p = [np.array(q["p"]) for q in prof["064"]["pts"] if q["ok"] and simulate(np.array(q["p"]), none)[0] < H][0]
g = np.round(np.arange(0.30, 1.19, 0.02), 2); sp = np.round(np.arange(0.30, 1.21, 0.01), 2)
fam = [onoff(lo, hi) for lo in g for hi in g if hi > lo + 0.01] + [integral(t, k) for k in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in sp] + \
      [proportional(t, k) for k in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in sp]
F = np.array([(mb, d) for ttp, d, mb in (simulate(p, r) for r in fam) if ttp >= H])
def tr(lo, hi, cnt):
    def rule(t, x, s):
        xp = s.get("x", x); s["x"] = x; on = s.get("on", True)
        if on and x <= lo: on = False
        elif not on and x >= hi: on = True; cnt["level"] += 1
        elif not on and x > lo and x - xp > 0.03: on = True; cnt["trend"] += 1
        s["on"] = on; return float(on)
    return rule
print("1.20 on the grid repr:", repr(float(np.round(np.arange(0.30, 1.51, 0.02), 2)[45])), "line", repr(v3.LINE))
for hi in (1.18, 1.19, 1.199, 1.20, 1.201, 1.21, 1.30, 1.50):
    best = None
    for lo in g[g < hi - 0.01]:
        cnt = {"level": 0, "trend": 0}; ttp, d, mb = simulate(p, tr(lo, hi, cnt))
        if ttp < H: continue
        ok = F[F[:, 0] <= mb + 1e-9]; sc = ok[:, 1].min() / d if len(ok) else np.inf
        if best is None or sc > best[0]: best = (sc, lo, dict(cnt), round(mb, 4))
    print(f"hi {hi:<6} best lo {best[1]:.2f} score {best[0]:.4f} switch-ons {best[2]} mean burden {best[3]}")
print("lo 0.70 across hi (ttp, score, switch-ons, x at the level switch-on):")
for hi in (1.18, 1.19, 1.195, 1.199, 1.1999, 1.20):
    cnt = {"level": 0, "trend": 0}; xs = []; r0 = tr(0.70, hi, cnt)
    def r(t, x, s, r0=r0):
        n = cnt["level"]; u = r0(t, x, s)
        if cnt["level"] > n: xs.append(round(x, 5))
        return u
    ttp, d, mb = simulate(p, r); ok = F[F[:, 0] <= mb + 1e-9]
    print(f"  hi {hi}: ttp {ttp} score {ok[:, 1].min() / d if ttp >= H and len(ok) else None} {cnt} level-on at x {xs}")
