"""Dose vs mean burden for one patient's best band fit: modulation (bench v2 family) vs on/off containment."""
import json, os, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import v3; from load import patients
pid = sys.argv[1] if len(sys.argv) > 1 else '015'
prof = {r['pid']: r for r in json.load(open(os.path.join(HERE, '..', 'profile_all.json')))}[pid]
q = min((q for q in prof['pts'] if q['ok']), key=lambda q: q['cost'])
mods = ['mod_%.2f' % x for x in v3.TARGETS]
hys = ['hys_%.2f_%.2f' % (a, b) for a in np.arange(0.1, 1.01, 0.1) for b in np.arange(0.2, 1.16, 0.1) if b > a + 0.05]
rules = ['mtd', 'adaptive50', 'cycle1'] + mods + hys
ttp, dose, mb = v3.run(np.array([q['p']]), rules, v3.cycle1(patients()[pid])); ttp, dose, mb = ttp[0], dose[0], mb[0]
fig, ax = plt.subplots(figsize=(7, 4.5), dpi=130)
for idx, lab, c, mk in ((range(3, 3 + len(mods)), 'weekly modulation (v2 frontier family)', '#2a6fdb', 'o'),
                        (range(3 + len(mods), len(rules)), 'on/off with two thresholds', '#e07b00', 's')):
    s = [i for i in idx if ttp[i] >= v3.H]; f = [i for i in idx if ttp[i] < v3.H]
    ax.scatter(mb[s], dose[s], c=c, marker=mk, s=22, label=lab + ' (survives 5 y)')
    ax.scatter(mb[f], dose[f], facecolors='none', edgecolors=c, marker=mk, s=22, alpha=.5, label=lab + ' (progresses)')
for i, lab in ((0, 'MTD'), (1, 'adaptive 50%'), (2, "replay patient's cycle 1")):
    ax.scatter(mb[i], dose[i], c='k', marker='*' if ttp[i] >= v3.H else 'x', s=90); ax.annotate(lab, (mb[i], dose[i]), xytext=(6, 4), textcoords='offset points', fontsize=8)
ax.set_xlabel('mean burden over the run (x starting burden)'); ax.set_ylabel('drug-days used in 5 years')
ax.set_title('Patient %s, best fit inside the band: drug used vs burden held' % pid, fontsize=10)
ax.legend(fontsize=7, loc='upper right'); fig.text(.01, .01, 'errata (AI agent), model fitted to Bruchovsky 2006 IAS data; in-silico, not medical advice', fontsize=6, color='#666')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'frontier_%s.png' % pid))
