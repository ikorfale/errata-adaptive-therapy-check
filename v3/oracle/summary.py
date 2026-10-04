"""Summary of the oracle runs (term = end-state constraint on, noterm = off). Later files override earlier ones per
(patient, fit, q). Prints per patient the max and median oracle_ratio over live fits and burden levels, and counts."""
import re, collections, statistics as st
R = {}
for term, files in ((1, ("oracle_term.out", "oracle_term2.out", "oracle_term3.out")), (0, ("oracle_noterm.out", "oracle_noterm2.out", "oracle_noterm3.out"))):
    for f in files:
        for l in open(f):
            m = re.match(r"(\d+) fit (\d+) (q\d+) B ([\d.]+) frontier ([\d.]+) oracle (\S+)(?: ratio ([\d.]+))?", l)
            if m: R[term, m[1], int(m[2]), m[3]] = None if m[6] == "None" else float(m[7])
            elif "no feasible" in l: R[term, l.split()[0], int(l.split()[2]), "none"] = "nofeas"
for term in (1, 0):
    print("=== end-state constraint", "ON" if term else "OFF")
    pids = sorted({k[1] for k in R if k[0] == term}); allr = []
    for pid in pids:
        rs = [v for k, v in R.items() if k[0] == term and k[1] == pid and isinstance(v, float)]
        unev = sum(1 for k, v in R.items() if k[0] == term and k[1] == pid and v is None)
        nof = sum(1 for k, v in R.items() if k[0] == term and k[1] == pid and v == "nofeas")
        allr += rs
        print(f"{pid}: n {len(rs)} max {max(rs):.4f} median {st.median(rs):.4f} >1.01: {sum(r > 1.01 for r in rs)} | start on the line (not evaluated) {unev} | fits with no feasible family rule {nof}")
    print(f"ALL: patients {len(pids)}, (fit, q) cells {len(allr)}, ratio > 1.01 in {sum(r > 1.01 for r in allr)}, > 1.10 in {sum(r > 1.10 for r in allr)}, max {max(allr):.4f}, median {st.median(allr):.4f}")
