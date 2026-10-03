"""claude-sonnet-scout 71618: relax the eligibility filter from >= 2 live fits to >= 1. Which patients join, and do the
15 fixed entries from pd_entry.py (not tuned here) win robustly on them against the finer frontier of maxmin064.py?"""
import json, os, numpy as np
from cm_cancer_q01_v3 import DATA, simulate, onoff, H, WIN
from cm_cancer_q01_v31 import integral, proportional
from pd_entry import ENTRIES
prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
none = lambda t, x, s: 0.0
grid = np.round(np.arange(0.30, 1.19, 0.02), 2); sp = np.round(np.arange(0.30, 1.21, 0.01), 2)
fam = [onoff(lo, hi) for lo in grid for hi in grid if hi > lo + 0.01] + \
      [integral(t, g) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in sp] + \
      [proportional(t, g) for g in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in sp]
nlive = {}
for pid, r in prof.items():
    fits = [np.array(q["p"]) for q in r["pts"] if q["ok"]]
    nlive[pid] = [p for p in fits if simulate(p, none)[0] < H]
print("live-fit counts over", len(prof), "patients:", {c: sum(1 for v in nlive.values() if len(v) == c) for c in sorted({len(v) for v in nlive.values()})})
wins = 0
for pid in sorted(p for p, v in nlive.items() if len(v) == 1):
    p = nlive[pid][0]
    F = np.array([(mb, d) for ttp, d, mb in (simulate(p, r) for r in fam) if ttp >= H])
    row = []
    for name, (f, a) in ENTRIES.items():
        ttp, d, mb = simulate(p, f(*a))
        if ttp < H: row.append((name, "FAIL")); continue
        ok = F[F[:, 0] <= mb + 1e-9] if len(F) else F
        row.append((name, float("inf") if not len(ok) else round(float(ok[:, 1].min() / d), 4)))
    sc = [s for _, s in row if s != "FAIL"]
    w = [n for n, s in row if s != "FAIL" and s > WIN]; wins += len(w)
    print(pid, "contained family rules", len(F), "| entries failing", sum(1 for _, s in row if s == "FAIL"),
          "| best entry score", max(sc) if sc else None, "| wins", w, flush=True)
print("robust wins on the one-live-fit patients:", wins)
