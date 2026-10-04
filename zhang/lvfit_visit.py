"""zenith 72077: is the day-0-threshold model's remaining error the visit grid? Same fits as lvfit.py (THR=base), but keep
the trajectory and score it on the observed process: (a) model crossing day rounded UP to the patient's next actual visit
after the stop (the day the clinic would have restarted), (b) rounded up to a fixed 28-day grid from the stop (no use of
future visit dates), (c) model PSA at the real restart visit vs the PSA recorded there. Rows go to lvfit_visit.json."""
import csv, json, numpy as np
from collections import defaultdict
import lvfit_core as C
P = defaultdict(list)
for r in csv.DictReader(open('../data/zhang/zhang_long.csv')):
    if r['arm'] == 'Adaptive': P[r['pid']].append((float(r['day']), float(r['psa']), int(float(r['abi']))))
rows = []
for pid, v in sorted(P.items()):
    v.sort(); d = np.array([x[0] for x in v]); p = np.array([x[1] for x in v]); a = np.array([x[2] for x in v])
    ch = [i for i in range(1, len(a)) if a[i] != a[i - 1]]
    stops = [i for i in ch if a[i] == 0]; starts = [i for i in ch if a[i] == 1]
    offs = [(s, next(t for t in starts if t > s)) for s in stops if any(t > s for t in starts)]
    for j in range(1, len(offs)):
        s, t = offs[j]; thr = float(p[0])
        dose = np.zeros(int(d[s]) + C.MAXOFF + 1)
        for i in range(len(d)):
            if d[i] <= d[s]: dose[int(d[i]):] = a[i]
        dose[int(d[s]):] = 0.0
        r = C.fit(d[:s + 1], p[:s + 1], dose)
        traj = C.sim(r.x, dose, len(dose))[int(d[s]):]
        hit = np.nonzero(traj >= thr)[0]; off_hat = int(hit[0]) if len(hit) else C.MAXOFF
        later = [x - d[s] for x in d[s + 1:] if x - d[s] >= off_hat]
        off_visit = float(later[0]) if later else float(off_hat)      # beyond the last visit: leave unrounded
        off_grid = float(28 * np.ceil(off_hat / 28)) if off_hat > 0 else 28.0
        off = float(d[t] - d[s]); last = float(d[offs[j - 1][1]] - d[offs[j - 1][0]])
        prev_visit_off = float(d[t - 1] - d[s])                         # last off-visit before restart
        rows.append(dict(pid=pid, j=j, off=off, model=off_hat, model_visit=off_visit, model_grid=off_grid, last=last,
                         thr=thr, psa_restart=float(p[t]), model_psa_at_restart=float(traj[int(off)]),
                         psa_prev_visit=float(p[t - 1]), prev_visit_off=prev_visit_off, beyond_last_visit=not later))
        print(json.dumps(rows[-1]), flush=True)
json.dump(rows, open('lvfit_visit.json', 'w'), indent=0)
