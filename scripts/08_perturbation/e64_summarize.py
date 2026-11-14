# -*- coding: utf-8 -*-
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

con = pd.read_csv(r"Enrichr_LINCS_L1000_Chem_Pert_Consensus_Sigs.csv")
down = con[con.term.str.endswith("Down")].copy()
up = con[con.term.str.endswith("Up")].copy()
sig_down = down[down.adjusted_pvalue < 0.05].sort_values("adjusted_pvalue")
print("significant reversal (Down) consensus terms:", len(sig_down))
top = sig_down.head(15).copy()
top["drug"] = top.term.str.replace(" Down", "")
print(top[["drug", "pvalue", "adjusted_pvalue", "combined_score", "overlap_genes"]].to_string(index=False))
top.to_csv("L1000_逆转型药物Top15.csv", index=False, encoding="utf-8-sig")

# 图：top 15 逆转药物 -log10 adj p
fig, ax = plt.subplots(figsize=(8, 6))
d = top.iloc[::-1]
colors = ["#c0392b" if d.drug.iloc[i].lower() == "luteolin" else "#3b6ea5"
          for i in range(len(d))]
ax.barh(d.drug, -np.log10(d.adjusted_pvalue), color=colors)
ax.set_xlabel("-log10(adjusted p-value)")
ax.set_title("L1000 consensus: drugs that reverse the inflammation-up signature")
plt.tight_layout(); plt.savefig("fig_L1000逆转药物.png", dpi=150)

# 专门确认 luteolin
print("\n=== Luteolin consensus ===")
print(con[con.term.str.lower().str.startswith("luteolin")].to_string(index=False))

# CRISPR KO 交叉
ko = pd.read_csv(r"Enrichr_LINCS_L1000_CRISPR_KO_Consensus_Sigs.csv")
print("\n=== CRISPR KO Down (top 8) ===")
print(ko[ko.term.str.endswith("Down")].head(8)[["term","adjusted_pvalue","overlap_genes"]].to_string(index=False))
