"""Step 4a: on which patients does the model beat 'replay cycle 1'? Uses forecast_all.json (step 3).
Features: how bad the baseline is (cycle 2 unlike cycle 1), training fit, hindsight resistant fraction,
number of held-out points. Rank-biserial style: compare medians for wins vs losses + Mann-Whitney."""
import json, numpy as np
from scipy.stats import mannwhitneyu, spearmanr
R = json.load(open('forecast_all.json'))
win = np.array([r['model'] < r['base'] for r in R])
feat = {
 'baseline rmse': [r['base'] for r in R],
 'train rmse': [r['train'] for r in R],
 'hindsight rmse': [r['hind'] for r in R],
 'R0/(S0+R0) hindsight': [r['hind_p']['R0'] / (r['hind_p']['S0'] + r['hind_p']['R0']) for r in R],
 'rR/rS hindsight': [r['hind_p']['rR'] / r['hind_p']['rS'] for r in R],
 'held-out points': [r['n_test'] for r in R],
 'all points': [r['n'] for r in R],
 'near-best fits': [r['n_good'] for r in R],
 'fit spread end (log hi/lo)': [np.log((r['end_hi'] + .1) / (r['end_lo'] + .1)) for r in R],
}
print('wins %d of %d' % (win.sum(), len(R)))
print('%-28s %10s %10s %8s' % ('feature', 'med win', 'med loss', 'MW p'))
for k, v in feat.items():
    v = np.array(v, float)
    p = mannwhitneyu(v[win], v[~win]).pvalue
    print('%-28s %10.3g %10.3g %8.3f' % (k, np.median(v[win]), np.median(v[~win]), p))
d = np.array([r['base'] - r['model'] for r in R])
for k in ('baseline rmse', 'R0/(S0+R0) hindsight', 'fit spread end (log hi/lo)'):
    s = spearmanr(feat[k], d); print('spearman(%s, base-model) %.2f p=%.3f' % (k, s.statistic, s.pvalue))
b = np.array(feat['baseline rmse']); q = np.quantile(b, [1/3, 2/3])
for lo, hi, nm in ((0, q[0], 'baseline good'), (q[0], q[1], 'middle'), (q[1], 9, 'baseline bad')):
    m = (b >= lo) & (b < hi); print('%-14s n=%2d model wins %d (%.0f%%)' % (nm, m.sum(), win[m].sum(), 100 * win[m].mean()))
