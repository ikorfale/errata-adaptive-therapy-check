# Paired comparison on the same 36 off-periods: fitted model (three restart levels) vs 'same as last' vs 'last x LOPO ratio'.
import re, math, statistics as st
def load(f): return [(l.split()[0], float(a), float(b), float(c), float(t)) for l in open(f) if ' cycle ' in l
                     for a, b, c, t in re.findall(r'off (\d+) model (\d+) last (\d+) thr ([\d.]+)', l)]
B, M, N = load('lvfit_base.out'), load('lvfit.out'), load('lvfit_min.out')
lr = {}
for pid, o, m, l, t in B: lr.setdefault(pid, []).append(math.log(o / l))
shr = [l * math.exp(st.median([x for q, v in lr.items() if q != pid for x in v])) for pid, o, m, l, t in B]
obs = [o for _, o, _, _, _ in B]
for name, pred in (('model, restart = day-0 PSA', [m for _, _, m, _, _ in B]), ('model, restart = mean earlier restart PSA', [m for _, _, m, _, _ in M]),
                   ('model, restart = min earlier restart PSA', [m for _, _, m, _, _ in N]), ('same as last', [l for _, _, _, l, _ in B]), ('last x LOPO ratio', shr)):
    e = [p - o for p, o in zip(pred, obs)]
    print(f"{name:45s} MAE {st.mean(map(abs, e)):6.1f} d | median |err| {st.median(map(abs, e)):5.0f} | signed median {st.median(e):+5.0f} | late {sum(x > 0 for x in e)}/36")
base_m = [m for _, _, m, _, _ in B]
print("model(day-0) closer than last x ratio in", sum(abs(m - o) < abs(s - o) for m, s, o in zip(base_m, shr, obs)), "/ 36;",
      "ties", sum(abs(m - o) == abs(s - o) for m, s, o in zip(base_m, shr, obs)))
print("restart level: day-0 PSA median %.2f vs mean earlier restart PSA median %.2f" % (st.median(t for *_, t in B), st.median(t for *_, t in M)))
