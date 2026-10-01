import json, numpy as np
R = json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'v3_all.json'))); H = 1825; n = len(R)
print('patients', n, 'band sets total', sum(r['n_band'] for r in R))
none_all = sum(all(t >= H for t in r['none_ttp']) for r in R); none_any = sum(any(t >= H for t in r['none_ttp']) for r in R)
print('no drug never crosses the line within 5 y: all band sets %d, some set %d' % (none_all, none_any))
print('MTD survives 5 y on all sets %d, on some %d' % (sum(all(t >= H for t in r['mtd_ttp']) for r in R), sum(any(t >= H for t in r['mtd_ttp']) for r in R)))
print('v2 saturation (every modulation target survives on every set): %d' % sum(r['every_mod_survives'] for r in R))
ok = [r for r in R if not any(t >= H for t in r['none_ttp'])]
print('usable (untreated tumour crosses the line on every set): %d' % len(ok))
g = [x for r in ok for x in r['frontier_gain']]
print('frontier: modulation-only dose / (modulation+on-off) dose at the same burden: median %.3f, >1.01 in %d/%d points; patients with any >1.01: %d/%d'
      % (np.median(g), sum(x > 1.01 for x in g), len(g), sum(any(x > 1.01 for x in r['frontier_gain']) for r in ok), len(ok)))
for rule in ('mtd', 'adaptive50', 'cycle1', 'mod_1.00', 'mod_0.50'):
    beat = fail = below = 0; mins = []
    for r in ok:
        sc = r['rules'][rule]['scores']
        if any(s[0] == 'FAIL' for s in sc): fail += 1; continue
        if any(s[0] == 'BELOW' for s in sc): below += 1; continue
        m = min(s[1] for s in sc); mins.append(m); beat += m > 1.01
    print('%-10s robust beat (>1.01 on every set) %2d | fails on some set %2d | below frontier range %2d | survivors median min-score %s'
          % (rule, beat, fail, below, '%.2f' % np.median(mins) if mins else '-'))
