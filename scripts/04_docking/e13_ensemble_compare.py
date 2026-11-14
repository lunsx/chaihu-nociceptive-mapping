# -*- coding: utf-8 -*-
import os, csv
from scipy.stats import spearmanr, kendalltau

# 5KIR（原 PTGS2 结构）
k5 = {"luteolin": -7.011, "quercetin": -6.603, "cyperotundone": -6.006,
      "kaempferol": -5.313, "paeoniflorin": -5.066, "alpha_cyperone": -5.541}
# 6COX（另一 PTGS2 结构）
k6 = {"alpha_cyperone": -8.901, "luteolin": -8.497, "quercetin": -8.068,
      "cyperotundone": -7.649, "kaempferol": -6.597, "paeoniflorin": -3.149}

ligs = list(k5.keys())
def rank(d):
    order = sorted(d, key=lambda x: d[x])  # 越负越靠前
    return {g: i + 1 for i, g in enumerate(order)}
r5, r6 = rank(k5), rank(k6)

rho = spearmanr([r5[g] for g in ligs], [r6[g] for g in ligs])
tau = kendalltau([r5[g] for g in ligs], [r6[g] for g in ligs])
print("Spearman rho: %.3f (p=%.4f)" % (rho.statistic, rho.pvalue))
print("Kendall tau: %.3f" % tau.statistic)

work = r"C:\dockwork\ensemble"
rows = []
print("\nligand            5KIR   rank | 6COX   rank | rank shift")
for g in sorted(ligs, key=lambda x: r6[x]):
    shift = r5[g] - r6[g]
    print(f"{g:16s} {k5[g]:6.2f} {r5[g]:4d} | {k6[g]:6.2f} {r6[g]:4d} | {shift:+d}")
    rows.append([g, k5[g], r5[g], k6[g], r6[g], shift])

with open(os.path.join(work, "ensemble_compare.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["ligand", "5KIR_affinity", "5KIR_rank", "6COX_affinity", "6COX_rank", "rank_shift"])
    w.writerows(rows)
    w.writerow([])
    w.writerow(["spearman_rho", round(rho.statistic, 3), "p", round(rho.pvalue, 4)])
    w.writerow(["kendall_tau", round(tau.statistic, 3)])
print("\nsaved ensemble_compare.csv")
