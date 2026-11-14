"""Add multiplicity-adjusted SMR p-values (Bonferroni, BH) to the SMR/HEIDI table."""
import csv
import os

WORKDIR = os.environ.get("CHAIHU_MR_WORKDIR",
                         os.path.dirname(os.path.abspath(__file__)))
path = os.path.join(WORKDIR, "smr_heidi_table.csv")

with open(path, newline="", encoding="utf-8") as fh:
    rows = list(csv.DictReader(fh))
fields = list(rows[0].keys())

p = [float(r["smr_p"]) for r in rows]
n = len(p)

# Bonferroni
bonf = [min(1.0, x * n) for x in p]

# Benjamini-Hochberg
order = sorted(range(n), key=lambda i: p[i])
bh = [0.0] * n
prev = 1.0
for rank in range(n, 0, -1):
    i = order[rank - 1]
    val = min(prev, p[i] * n / rank)
    bh[i] = val
    prev = val

for i, r in enumerate(rows):
    r["smr_p_bonf"] = round(bonf[i], 4)
    r["smr_p_bh"] = round(bh[i], 4)

for extra in ("smr_p_bonf", "smr_p_bh"):
    if extra not in fields:
        fields.append(extra)

with open(path, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f"n tests={n}, Bonferroni threshold={0.05 / n:.4f}")
print("nominal p<0.05:", [r["gene"] + "/" + r["outcome"] for r in rows
                          if float(r["smr_p"]) < 0.05])
print("Bonferroni p<0.05:", [r["gene"] + "/" + r["outcome"] for r in rows
                             if float(r["smr_p_bonf"]) < 0.05])
print("BH p<0.05:", [r["gene"] + "/" + r["outcome"] for r in rows
                     if float(r["smr_p_bh"]) < 0.05])
