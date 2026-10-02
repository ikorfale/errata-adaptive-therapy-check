"""Holdout cohort check (errata, 2026-10-02): how many Shaw et al. patients would the bench be able to rank?
Same v3 scorer (v3.score_patient) on the band sets of profile_shaw.json. A set is live if the untreated model
tumour crosses the progression line (1.2 N0) within H; a patient is eligible with >= 2 live sets (ruling 10-01).
Also prints the sha256 of the band file, so the holdout can be committed before any rule is scored on it."""
import glob, hashlib, json, os, sys
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); sys.path.insert(0, os.path.dirname(here))
import v3, load
D = os.path.join(os.path.dirname(here), 'data/bruchovsky/dataTanaka/Shaw_et_al')
PF = os.path.join(os.path.dirname(here), 'profile_shaw.json')
rows = {os.path.basename(p)[7:-4]: load.load(p) for p in sorted(glob.glob(D + '/patient*.txt'))}
prof = json.load(open(PF))
print('band file sha256', hashlib.sha256(open(PF, 'rb').read()).hexdigest(), '| patients fitted', len(prof), 'of', len(rows))
elig = 0; dist = {}
for r in prof:
    s = v3.score_patient(r['pid'], rows[r['pid']], r)
    live = sum(t < v3.H for t in s['none_ttp']); n = s['n_band']; dist[live] = dist.get(live, 0) + 1
    elig += live >= 2
    print(r['pid'], f'band {n}, live {live}, untreated ttp {s["none_ttp"]}')
print('live sets per patient:', dict(sorted(dist.items())), '| eligible (>=2 live):', elig, 'of', len(prof))
