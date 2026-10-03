# Zhang adaptive arm: split each patient's visits into treatment cycles (on-run then off-run of abi),
# then ask how well "the next cycle is like the last one" predicts the next cycle's on- and off-durations.
import csv, numpy as np
from collections import defaultdict
P = defaultdict(list)
for r in csv.DictReader(open('../data/zhang/zhang_long.csv')):
    if r['arm'] == 'Adaptive': P[r['pid']].append((float(r['day']), float(r['psa']), int(float(r['abi']))))
cyc = {}
for pid, v in sorted(P.items()):
    v.sort(); runs = []                      # (state, start_day) at each change of abi
    for d, psa, a in v:
        if not runs or runs[-1][0] != a: runs.append((a, d))
    end = v[-1][0]
    spans = [(a, (runs[i + 1][1] if i + 1 < len(runs) else end) - s, i + 1 < len(runs)) for i, (a, s) in enumerate(runs)]
    c = []                                   # complete cycles: an on-span then an off-span, both closed
    for i in range(len(spans) - 1):
        if spans[i][0] == 1 and spans[i + 1][0] == 0 and spans[i][2] and spans[i + 1][2]:
            c.append((spans[i][1], spans[i + 1][1]))
    cyc[pid] = c
    print(pid, 'visits', len(v), 'days', int(end), 'complete cycles', len(c), [(int(a), int(b)) for a, b in c])
for k, name in ((0, 'on'), (1, 'off')):
    err_last, err_mean, n = [], [], 0
    for pid, c in cyc.items():
        for j in range(1, len(c)):
            y = c[j][k]; err_last.append(abs(c[j - 1][k] - y)); err_mean.append(abs(np.mean([x[k] for x in c[:j]]) - y))
    print('%s-duration: n %d next cycles | MAE last-cycle %.1f d | MAE mean-of-previous %.1f d | median %s %.0f d'
          % (name, len(err_last), np.mean(err_last), np.mean(err_mean), name, np.median([x[k] for c in cyc.values() for x in c])))
# Do off-periods shrink? log-ratio of successive off-durations; leave-one-patient-out shrink rule.
lr = {pid: [np.log(c[j][1] / c[j - 1][1]) for j in range(1, len(c))] for pid, c in cyc.items()}
allr = [x for v in lr.values() for x in v]
print('off log-ratio next/last: n %d, median ratio %.2f, shrinks in %d/%d' % (len(allr), np.exp(np.median(allr)), sum(x < 0 for x in allr), len(allr)))
e_shr = []
for pid, c in cyc.items():
    other = [x for q, v in lr.items() if q != pid for x in v]; f = np.exp(np.median(other))
    e_shr += [abs(c[j - 1][1] * f - c[j][1]) for j in range(1, len(c))]
print('off-duration MAE, last-cycle x LOPO median ratio: %.1f d (n %d)' % (np.mean(e_shr), len(e_shr)))
