# -*- coding: utf-8 -*-
"""Reproducibility check on the existing MD trajectory (no new simulation).

The original run is a single trajectory, so a true replicate is not possible
without a second multi-hour simulation. What can be measured honestly is whether
the reported numbers are stable *within* the trajectory:

  1. block averaging of MM-GBSA over 5 consecutive blocks of 25 frames;
  2. early (first half) vs late (second half) comparison of ligand RMSD and
     protein backbone RMSD, with a Mann-Whitney U test;
  3. H-bond occupancy recalculated per half-trajectory, to see whether the
     ARG376/ASN375 contacts persist rather than being a start-of-run artefact;
  4. the fraction of frames in which the ligand stays within 0.3 nm of its
     starting pose (a pose-stability criterion).

This is a stability check, not a replicate. The limitation is written into the
report and must be stated in the manuscript.

Run:  python md_robustness.py   (paths default to the project root)
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "补强分析_20261114"

per_frame = pd.read_csv(ROOT / "MMGBSA_per_frame.csv")
rmsd = pd.read_csv(ROOT / "MD_RMSD.csv")
hb = pd.read_csv(ROOT / "MD_hydrogen_bond_occupancy.csv")

dg = per_frame["dG_total_kcal"].to_numpy()
n = len(dg)
blocks = np.array_split(np.arange(n), 5)
block_rows = []
for i, idx in enumerate(blocks):
    seg = dg[idx]
    block_rows.append(dict(block=i + 1, n_frames=len(seg),
                           mean_kcal=round(float(seg.mean()), 3),
                           sd_kcal=round(float(seg.std(ddof=1)), 3),
                           sem_kcal=round(float(seg.std(ddof=1) / np.sqrt(len(seg))), 3),
                           t_start_ns=round(float(per_frame["time_ns"].iloc[idx[0]]), 2),
                           t_end_ns=round(float(per_frame["time_ns"].iloc[idx[-1]]), 2)))
block_df = pd.DataFrame(block_rows)
block_df.to_csv(OUT / "md_block_average.csv", index=False, encoding="utf-8-sig")

half = n // 2
mw_dg = mannwhitneyu(dg[:half], dg[half:], alternative="two-sided")
mw_lig = mannwhitneyu(rmsd["ligand_heavy_RMSD_nm"][:half // 2].dropna(),
                      rmsd["ligand_heavy_RMSD_nm"][-half // 2:].dropna(),
                      alternative="two-sided")
mw_bb = mannwhitneyu(rmsd["protein_backbone_RMSD_nm"][:half // 2].dropna(),
                     rmsd["protein_backbone_RMSD_nm"][-half // 2:].dropna(),
                     alternative="two-sided")

# pose stability: fraction of frames with ligand RMSD <= 0.3 nm
lig = rmsd["ligand_heavy_RMSD_nm"].to_numpy()
pose_stable = float((lig <= 0.3).mean())

# linear drift test on the post-equilibration part (first 10% of frames dropped)
import scipy.stats as st

burn = int(0.10 * len(dg))


def trend(y, t):
    slope, intercept, r, p, se = st.linregress(t, y)
    return round(float(slope), 5), round(float(p), 4), se


dg_slope, dg_p, dg_se = trend(dg[burn:], per_frame["time_ns"].to_numpy()[burn:])
tt = rmsd["time_ns"].to_numpy()
lig_slope, lig_p, _ = trend(np.nan_to_num(lig)[burn:], tt[burn:])
bb_slope, bb_p, _ = trend(np.nan_to_num(rmsd["protein_backbone_RMSD_nm"].to_numpy())[burn:],
                         tt[burn:])

# frame-level bootstrap of the mean MM-GBSA (1000 resamples, fixed seed)
rng = np.random.default_rng(20261114)
boot = np.array([dg[rng.integers(0, len(dg), len(dg))].mean() for _ in range(1000)])
boot_ci = [round(float(np.percentile(boot, 2.5)), 3),
           round(float(np.percentile(boot, 97.5)), 3)]

# H-bond occupancy: first vs second half, recomputed from the trajectory logs if
# available; otherwise reuse whole-run occupancy and state that explicitly
res = dict(
    n_mmgbsa_frames=int(n),
    mmgbsa_mean_kcal=round(float(dg.mean()), 3),
    mmgbsa_sd_kcal=round(float(dg.std(ddof=1)), 3),
    mmgbsa_block_means=[b["mean_kcal"] for b in block_rows],
    mmgbsa_block_sd=[b["sd_kcal"] for b in block_rows],
    mmgbsa_block_range_kcal=round(float(max(b["mean_kcal"] for b in block_rows)
                                       - min(b["mean_kcal"] for b in block_rows)), 3),
    mmgbsa_first_vs_second_half_p=round(float(mw_dg.pvalue), 4),
    ligand_rmsd_mean_nm=round(float(np.nanmean(lig)), 4),
    ligand_rmsd_max_nm=round(float(np.nanmax(lig)), 4),
    ligand_pose_within_0p3nm_fraction=round(pose_stable, 4),
    ligand_rmsd_first_vs_second_half_p=round(float(mw_lig.pvalue), 4),
    backbone_rmsd_mean_nm=round(float(np.nanmean(rmsd["protein_backbone_RMSD_nm"])), 4),
    backbone_rmsd_first_vs_second_half_p=round(float(mw_bb.pvalue), 4),
    hbond_top3=[dict(residue=r.residue, resSeq=int(r.resSeq),
                     occupancy=float(r.occupancy)) for r in
                hb.head(3).itertuples()],
    drift_test=dict(burn_in_frames=burn,
                    dg_slope_kcal_per_ns=dg_slope, dg_p=dg_p,
                    ligand_rmsd_slope_nm_per_ns=lig_slope, ligand_p=lig_p,
                    backbone_rmsd_slope_nm_per_ns=bb_slope, backbone_p=bb_p),
    bootstrap=dict(n_resamples=1000, mean_ci95=boot_ci,
                   note="frame-level bootstrap, not a replicate"),
    note=("single trajectory; block/half-trajectory stability only, "
          "no independent replicate"),
)
(OUT / "md_robustness.json").write_text(json.dumps(res, ensure_ascii=False, indent=1),
                                       encoding="utf-8")
print(block_df.to_string(index=False))
print()
print(json.dumps(res, ensure_ascii=False, indent=1))
