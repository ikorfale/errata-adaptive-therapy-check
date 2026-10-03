# TrialPatientData.xlsx (Zhang et al. 2022 eLife, CC0, cunninghamjj repo) -> long CSV: arm,pid,day,psa,abi
import openpyxl, csv, warnings; warnings.simplefilter('ignore')
wb = openpyxl.load_workbook('TrialPatientData.xlsx', read_only=True)
out = csv.writer(open('zhang_long.csv', 'w', newline='')); out.writerow(['arm', 'pid', 'day', 'psa', 'abi'])
for ws in wb.worksheets:
    rows = list(ws.iter_rows(values_only=True)); head = rows[0]
    for c, pid in enumerate(head):
        if not pid: continue
        n = 0
        for r in rows[2:]:
            if c + 2 >= len(r) or r[c] is None or r[c + 1] is None: continue
            out.writerow([ws.title, pid, r[c], r[c + 1], r[c + 2]]); n += 1
        print(ws.title, pid, n, 'visits')
