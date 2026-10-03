# claude-sonnet-scout (71561): separate "wins on every live fit" from "wins on the mean over live fits".
# A FAIL (tumour crosses the line on a fit) has no dose ratio; report two readings:
#   strict-mean: a FAIL anywhere still disqualifies (safety is not averaged); mean over OK scores
#   loose-mean : FAIL fits dropped, mean over the fits where the rule held (an optimistic upper bound)
import json, numpy as np
d = json.load(open('CM-CANCER-Q01-v31.json')); WIN = 1.01
P = {k: v for k, v in d['patients'].items() if v['eligible']}
rules = list(d['summary']['entries'])
print('eligible patients', len(P))
for r in rules:
    every = smean = lmean = 0; best = []
    for pid, v in P.items():
        live = [s[r] for s in v['sets'] if s[r][0] != 'VACUOUS']
        ok = [s[1] for s in live if s[0] == 'OK' and s[1] is not None]
        nfail = sum(s[0] == 'FAIL' for s in live); nbelow = sum(s[0] == 'OK' and s[1] is None for s in live)  # below frontier range
        if ok and len(ok) == len(live) and min(ok) > WIN: every += 1
        if ok and len(ok) == len(live) and np.mean(ok) > WIN: smean += 1
        if ok and np.mean(ok) > WIN: lmean += 1
        if ok: best.append((round(float(np.mean(ok)), 3), round(max(ok), 3), pid, len(ok), nfail, nbelow))
    best.sort(reverse=True)
    print('%-13s every-fit %d | mean, no fail %d | mean over held fits %d | top (mean,max,pid,ok,fail,below) %s'
          % (r, every, smean, lmean, best[:2]))
allok = [s[r][1] for v in P.values() for s in v['sets'] for r in rules if s[r][0] == 'OK' and s[r][1] is not None]
print('all OK live-fit scores: n %d, max %.3f, >1.01: %d' % (len(allok), max(allok), sum(x > WIN for x in allok)))
