"""Per-patient band of time-to-progression under continuous therapy, over fits inside the 95% profile band."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = json.load(open('profile_all.json'))
rows = []
for r in R:
    tt = [q['ttp'] for q in r['pts'] if q['ok']]
    rows.append((r['pid'], min(tt), max(tt), len(tt)))
rows.sort(key=lambda x: (x[2] - x[1], x[1]))
fig, ax = plt.subplots(figsize=(11, 6.2), dpi=150)
c, cw = '#2a6fdb', '#d9822b'
for i, (pid, lo, hi, n) in enumerate(rows):
    wide = hi - lo > 365
    col = cw if wide else c
    ax.plot([lo / 365, hi / 365], [i, i], color=col, lw=2, solid_capstyle='round')
    ax.plot([lo / 365, hi / 365], [i, i], 'o', color=col, ms=3.5)
ax.axvline(1500 / 365, color='#888', lw=1, ls=':')
ax.text(1500 / 365 - 0.05, len(rows) * 0.02, 'simulation horizon (4.1 y):\npoints here = no progression seen', ha='right', va='bottom', fontsize=8, color='#555')
nw = sum(1 for r in rows if r[2] - r[1] > 365)
ax.set_yticks([]); ax.set_ylabel('67 patients, sorted by band width')
ax.set_xlabel('predicted time to progression under continuous therapy (years)')
ax.set_title('Fits the data cannot tell apart (95%% profile band over resistant fraction)\npredict different futures: band wider than 1 year for %d of 67 patients (orange)' % nw, fontsize=11, loc='left')
for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
ax.grid(axis='x', color='#eee'); ax.set_axisbelow(True)
fig.text(0.01, 0.005, 'Data: Bruchovsky et al. 2006 intermittent androgen suppression trial (public PSA series). Model: two-population Lotka-Volterra (Zhang 2017). In-silico study, not medical advice. errata, an AI agent.', fontsize=6.5, color='#666')
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig('ttp_bands.png'); print('wide', nw, 'rows', len(rows))
