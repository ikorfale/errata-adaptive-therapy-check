"""Load the public Bruchovsky et al. (2006) IAS trial data (dataTanaka.zip, nicholasbruchovsky.com).
Columns: patient, date, CPA, LEU, PSA, testosterone, cycle, on(1)/off(0), day, day2."""
import csv, glob, os
D = os.path.join(os.path.dirname(__file__), 'data/bruchovsky/dataTanaka/Bruchovsky_et_al')
def num(x):
    try: return float(x.strip('"'))
    except ValueError: return None
def load(path):
    rows = []
    for r in csv.reader(open(path)):
        if len(r) < 9: continue
        rows.append(dict(day=num(r[8]), psa=num(r[4]), testo=num(r[5]), cycle=num(r[6]), on=num(r[7])))
    rows = sorted((r for r in rows if r['day'] is not None), key=lambda r: r['day'])  # 016, 037 are stored out of order
    d0 = rows[0]['day']
    for r in rows: r['t'] = r['day'] - d0
    return rows
def patients():
    return {os.path.basename(p)[7:10]: load(p) for p in sorted(glob.glob(D + '/patient*.txt'))}
if __name__ == '__main__':
    P = patients(); import statistics as s
    n_psa = [sum(r['psa'] is not None for r in v) for v in P.values()]
    cyc = [max(r['cycle'] or 0 for r in v) for v in P.values()]
    span = [v[-1]['t']/365 for v in P.values()]
    print('patients', len(P), 'PSA points median', s.median(n_psa), 'min', min(n_psa), 'max', max(n_psa))
    print('cycles: ', {c: cyc.count(c) for c in sorted(set(cyc))})
    print('follow-up years median %.1f max %.1f' % (s.median(span), max(span)))
    print('patients with >=2 cycles', sum(c >= 2 for c in cyc))
