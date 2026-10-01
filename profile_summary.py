"""Summarise profile_all.json (step 4b): how loosely do the training data pin the future?"""
import json, numpy as np
R = json.load(open('profile_all.json'))
F = {r['pid']: r for r in json.load(open('forecast_all.json'))}
print('patients', len(R))
nok = np.array([r['n_ok'] for r in R]); print('grid points inside 95%% band (of 9): median %d, >=3 for %d, all 9 for %d' % (np.median(nok), (nok >= 3).sum(), (nok == 9).sum()))
for key, lab in (('ok', '95% band'), ('ok20', 'cost within 20%')):
    sp, cens = [], 0
    for r in R:
        tt = [q['ttp'] for q in r['pts'] if q[key]]
        sp.append(max(tt) - min(tt)); cens += (max(tt) == 1500) and (min(tt) < 1500)
    sp = np.array(sp)
    print('%-16s TTP spread under continuous therapy (days): median %d, >180 d for %d, >365 d for %d; band spans "never in 4 y" and a finite TTP for %d' % (lab, np.median(sp), (sp > 180).sum(), (sp > 365).sum(), cens))
fcs = np.array([r['fc_hi'] - r['fc_lo'] for r in R]); print('forecast RMSE spread inside 95%% band: median %.2f, >0.3 for %d' % (np.median(fcs), (fcs > .3).sum()))
win_any = sum(r['fc_lo'] < F[r['pid']]['base'] for r in R if r['pid'] in F); win_all = sum(r['fc_hi'] < F[r['pid']]['base'] for r in R if r['pid'] in F)
print('beats replay baseline with the best in-band fit (chosen with hindsight, optimistic): %d; with every in-band fit: %d (of %d)' % (win_any, win_all, len(R)))
