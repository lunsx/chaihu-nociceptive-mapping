# -*- coding: utf-8 -*-
"""Mouse MASLD Visium: donor-level case-control descriptive table.

The 4a34293d dataset has one SD (normal) and two WD (MASLD) donors.
We report per-donor mean expression and positive-spot fraction for the core
genes, plus the WD-minus-SD contrast. No p-values: the normal arm has n=1.

Path configuration: set CHAIHU_SPATIAL_DIR to the directory holding
mouse_masld_visium.h5ad; the cloud analysis path is the default.
"""
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

_SPATIAL = Path(os.environ.get(
    "CHAIHU_SPATIAL_DIR", "/home/user/compute_work/spatial"))
F = str(_SPATIAL / "mouse_masld_visium.h5ad")
OUT = _SPATIAL
GENES = ["Il1b", "Tnf", "Ptgs2", "Il6", "Ptges", "Ptger3", "Calcr l",
         "Tacr1", "Trpv1", "Fos", "Spp1", "Hmox1", "Icam1", "Ccl2"]
GENES = ["Il1b", "Tnf", "Ptgs2", "Il6", "Ptges", "Ptger3", "Calcr",
         "Tacr1", "Trpv1", "Fos", "Spp1", "Hmox1", "Icam1", "Ccl2"]

a = ad.read_h5ad(F)
sym2ens = dict(zip(a.var["feature_name"].astype(str), a.var_names))
donors = a.obs["donor_id"].astype(str).to_numpy()
disease = a.obs["disease"].astype(str).to_numpy()
rows = []
for g in GENES:
    ens = sym2ens.get(g)
    if ens is None:
        continue
    x = a[:, ens].X
    x = np.asarray(x.toarray()).ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
    per_donor = {}
    for d in sorted(set(donors)):
        m = donors == d
        per_donor[d] = (float(x[m].mean()), float((x[m] > 0).mean()), int(m.sum()))
    sd = [d for d in per_donor if d.startswith("SD")]
    wd = [d for d in per_donor if d.startswith("WD")]
    sd_mean = np.mean([per_donor[d][0] for d in sd])
    wd_mean = np.mean([per_donor[d][0] for d in wd])
    sd_pct = np.mean([per_donor[d][1] for d in sd]) * 100
    wd_pct = np.mean([per_donor[d][1] for d in wd]) * 100
    rows.append({
        "gene": g,
        "SD_mean": round(sd_mean, 4),
        "WD_mean": round(wd_mean, 4),
        "WD_minus_SD": round(wd_mean - sd_mean, 4),
        "fold_change": (round((wd_mean + 1e-6) / (sd_mean + 1e-6), 2)
                        if sd_mean > 0 else "n/a (SD=0)"),
        "SD_pct_positive_spots": round(sd_pct, 2),
        "WD_pct_positive_spots": round(wd_pct, 2),
        "SD_donor_n": len(sd),
        "WD_donor_n": len(wd),
    })

out = pd.DataFrame(rows)
out.to_csv(OUT / "空间_小鼠病例对照_供体级.csv", index=False, encoding="utf-8-sig")
print(out.to_string(index=False))
print("\nNote: normal arm n=1 donor; descriptive contrast only, no inferential test.")
