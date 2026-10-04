"""Score lvfit_visit.json: pooled MAE, per-patient mean MAE, patient bootstrap of (model - last x LOPO ratio), and a test
of zenith's premise (restart = first visit after PSA crosses day-0 PSA)."""
import json, math, statistics as st, numpy as np
R = json.load(open('lvfit_visit.json')); obs = [r['off'] for r in R]
lr = {}
for r in R: lr.setdefault(r['pid'], []).append(math.log(r['off'] / r['last']))
shr = [r['last'] * math.exp(st.median([x for q, v in lr.items() if q != r['pid'] for x in v])) for r in R]
import csv
from collections import defaultdict
P = defaultdict(list)
for x in csv.DictReader(open('../data/zhang/zhang_long.csv')):
    if x['arm'] == 'Adaptive': P[x['pid']].append((float(x['day']), int(float(x['abi']))))
STOP = {}
for pid, v in P.items():
    v.sort(); d = [x[0] for x in v]; a = [x[1] for x in v]; ch = [i for i in range(1, len(a)) if a[i] != a[i - 1]]
    stops = [i for i in ch if a[i] == 0]; starts = [i for i in ch if a[i] == 1]
    offs = [(s, next(t for t in starts if t > s)) for s in stops if any(t > s for t in starts)]
    for j, (s, t) in enumerate(offs): STOP[pid, j] = (d[s], d[s + 1:])
def nextvisit(r, x):   # round a predicted off-day up to the patient's next actual visit after the stop (as model_visit)
    s, later = STOP[r['pid'], r['j']]; c = [v - s for v in later if v - s >= x]; return c[0] if c else x
assert all(nextvisit(r, r['model']) == r['model_visit'] for r in R)
preds = {'model, crossing day': [r['model'] for r in R], 'model, next actual visit': [r['model_visit'] for r in R],
         'model, 28-day grid': [r['model_grid'] for r in R], 'same as last': [r['last'] for r in R],
         'last x LOPO ratio': shr, 'last x LOPO ratio, 28-day grid': [28 * math.ceil(s / 28) for s in shr],
         'last x LOPO ratio, next actual visit': [nextvisit(r, s) for r, s in zip(R, shr)]}
pids = sorted(set(r['pid'] for r in R)); rng = np.random.default_rng(0)
for name, pr in preds.items():
    e = np.array([p - o for p, o in zip(pr, obs)])
    pp = [np.mean(abs(e[[i for i, r in enumerate(R) if r['pid'] == q]])) for q in pids]
    print(f"{name:32s} MAE {np.mean(abs(e)):6.1f} | per-patient mean MAE {np.mean(pp):6.1f} | median |err| {np.median(abs(e)):5.0f} | signed median {np.median(e):+5.0f}")
def boot(a, b, n=10000):
    ea, eb = np.abs(np.array(a) - obs), np.abs(np.array(b) - obs); idx = {q: [i for i, r in enumerate(R) if r['pid'] == q] for q in pids}
    out = []
    for _ in range(n):
        ii = [i for q in rng.choice(pids, len(pids)) for i in idx[q]]; out.append(ea[ii].mean() - eb[ii].mean())
    lo, hi = np.percentile(out, [2.5, 97.5]); return np.mean(ea) - np.mean(eb), lo, hi
for a, b in (('model, crossing day', 'last x LOPO ratio'), ('model, next actual visit', 'last x LOPO ratio, next actual visit'), ('model, 28-day grid', 'last x LOPO ratio, 28-day grid')):
    d, lo, hi = boot(preds[a], preds[b])
    print(f"MAE({a}) - MAE({b}) = {d:+.1f} d, patient-bootstrap 95% CI [{lo:+.1f}, {hi:+.1f}]")
print("cycles per patient:", {q: sum(r['pid'] == q for r in R) for q in pids})
# premise: restart PSA >= day-0 PSA, and the previous off-visit below it (crossing between visits)
above = sum(r['psa_restart'] >= r['thr'] for r in R); prev_below = sum(r['psa_prev_visit'] < r['thr'] for r in R)
both = sum(r['psa_restart'] >= r['thr'] and r['psa_prev_visit'] < r['thr'] for r in R)
print(f"restart PSA >= day-0 PSA: {above}/36 | previous visit below day-0 PSA: {prev_below}/36 | both (crossing between the two visits): {both}/36")
print("restart PSA / day-0 PSA: median %.2f, range %.2f-%.2f" % (st.median(r['psa_restart'] / r['thr'] for r in R), min(r['psa_restart'] / r['thr'] for r in R), max(r['psa_restart'] / r['thr'] for r in R)))
lg = [math.log(r['model_psa_at_restart'] / r['psa_restart']) for r in R]
print("model PSA at the real restart visit / observed restart PSA: median %.2f, within x2 in %d/36" % (math.exp(st.median(lg)), sum(abs(x) < math.log(2) for x in lg)))
print("restart at the first visit after the stop (no off-visit in between):", sum(r['prev_visit_off'] == 0 for r in R), "/36")
print("model crossing beyond the patient's last visit (left unrounded):", sum(r['beyond_last_visit'] for r in R))
