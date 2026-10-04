"""zenith-claude 72182: score forecasts against the PSA crossing itself, as an interval. For each of the 36 cycles the
crossing of day-0 PSA lies between the last off-visit below it and the first off-or-restart visit at or above it
(days after the stop). Error = 0 inside the bracket, else distance to the nearer end. Cycles with no visit at or above
day-0 PSA before restart have no bracket and are skipped (counted). Forecasts: model crossing day (lvfit_visit.json),
last x LOPO ratio fitted on restart dates (as posted), and last-bracket-midpoint x LOPO ratio fitted on bracket midpoints."""
import csv, json, math, statistics as st, numpy as np
from collections import defaultdict
R = json.load(open('lvfit_visit.json'))
P = defaultdict(list)
for x in csv.DictReader(open('../data/zhang/zhang_long.csv')):
    if x['arm'] == 'Adaptive': P[x['pid']].append((float(x['day']), float(x['psa']), int(float(x['abi']))))
BR = {}
for pid, v in P.items():
    v.sort(); d = [x[0] for x in v]; p = [x[1] for x in v]; a = [x[2] for x in v]; thr = p[0]
    ch = [i for i in range(1, len(a)) if a[i] != a[i - 1]]; stops = [i for i in ch if a[i] == 0]; starts = [i for i in ch if a[i] == 1]
    offs = [(s, next(t for t in starts if t > s)) for s in stops if any(t > s for t in starts)]
    for j, (s, t) in enumerate(offs):
        up = [i for i in range(s + 1, t + 1) if p[i] >= thr]
        if not up: BR[pid, j] = None; continue
        i = up[0]; BR[pid, j] = (d[i - 1] - d[s], d[i] - d[s])
def err(x, b): lo, hi = b; return 0.0 if lo <= x <= hi else min(abs(x - lo), abs(x - hi))
rows = [r for r in R if BR[r['pid'], r['j']] is not None]
print(f"cycles with a bracket: {len(rows)} of {len(R)}; bracket width median {st.median(BR[r['pid'], r['j']][1] - BR[r['pid'], r['j']][0] for r in rows):.0f} d")
lr = defaultdict(list)
for r in R: lr[r['pid']].append(math.log(r['off'] / r['last']))
shr = {(r['pid'], r['j']): r['last'] * math.exp(st.median([x for q, v in lr.items() if q != r['pid'] for x in v])) for r in R}
mid = lambda pid, j: None if BR.get((pid, j)) is None else sum(BR[pid, j]) / 2
lm = defaultdict(list)                       # ratio of successive bracket midpoints, for a naive rule trained on the crossing
for r in R:
    a, b = mid(r['pid'], r['j'] - 1), mid(r['pid'], r['j'])
    if a and b: lm[r['pid']].append(math.log(b / a))
F = {'model crossing day': lambda r: r['model'], 'last x LOPO ratio (restart-trained)': lambda r: shr[r['pid'], r['j']],
     'last bracket midpoint x LOPO ratio': lambda r: (mid(r['pid'], r['j'] - 1) * math.exp(st.median([x for q, v in lm.items() if q != r['pid'] for x in v]))
                                                      if mid(r['pid'], r['j'] - 1) else None)}
pids = sorted({r['pid'] for r in rows}); rng = np.random.default_rng(0); E = {}
for name, f in F.items():
    sub = [(r, f(r)) for r in rows if f(r) is not None]; e = [err(x, BR[r['pid'], r['j']]) for r, x in sub]
    E[name] = {(r['pid'], r['j']): v for (r, _), v in zip(sub, e)}
    print(f"{name:40s} n {len(e)} | mean interval error {st.mean(e):5.1f} d | median {st.median(e):4.0f} | inside bracket {sum(v == 0 for v in e)}/{len(e)}")
def boot(a, b):
    keys = sorted(set(E[a]) & set(E[b])); idx = defaultdict(list)
    for k in keys: idx[k[0]].append(k)
    ps = sorted(idx); d0 = st.mean(E[a][k] - E[b][k] for k in keys); out = []
    for _ in range(10000):
        ks = [k for q in rng.choice(ps, len(ps)) for k in idx[q]]; out.append(np.mean([E[a][k] - E[b][k] for k in ks]))
    lo, hi = np.percentile(out, [2.5, 97.5]); return len(keys), len(ps), d0, lo, hi
for b in ('last x LOPO ratio (restart-trained)', 'last bracket midpoint x LOPO ratio'):
    n, np_, d0, lo, hi = boot('model crossing day', b)
    print(f"model - {b}: {d0:+.1f} d on {n} cycles / {np_} patients, patient-bootstrap 95% CI [{lo:+.1f}, {hi:+.1f}]")
