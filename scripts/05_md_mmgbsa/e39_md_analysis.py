# -*- coding: utf-8 -*-
"""
e39  MD 轨迹分析（5KIR-luteolin）
产出：MD_RMSD.csv / MD_RMSF.csv / MD_hydrogen_bonds.csv / MD_lig_prot_distance.csv
      以及 Fig8_MD稳定性.png（投稿图集）
依赖：mdtraj；读取 C:\\dockwork\\md\\modeller.pkl 与 traj.dcd
说明：脚本可在 MD 跑完后直接运行；若提前运行，读到的是当前已有帧。
"""
import os, pickle
import numpy as np, pandas as pd
import mdtraj as md

MD=r"C:\dockwork\md"
ROOT=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
OUTDIR=os.path.join(ROOT,"投稿图集"); os.makedirs(OUTDIR,exist_ok=True)

m=pickle.load(open(os.path.join(MD,"modeller.pkl"),"rb"))
top=md.Topology.from_openmm(m["topology"])
lig_atoms=list(m["lig_idx"])
n_prot=int(m["n_prot"])
prot_atoms=list(range(0,n_prot))
backbone=top.select("backbone and chainid 0")
ca=top.select("name CA and chainid 0")
def is_heavy(i): return top.atom(i).element.symbol!="H"
lig_heavy=[i for i in lig_atoms if is_heavy(i)]

tr_all=md.load(os.path.join(MD,"traj.dcd"),top=top)
# DCD 前 6 帧（step 10000..60000）是 e38 平衡/约束阶段写入，生产帧从索引 6（step 70000）起
SKIP=6
t_all=np.arange(tr_all.n_frames)*0.04   # 40 ps/帧
prod_idx=list(range(SKIP,tr_all.n_frames))
tr=tr_all.slice(prod_idx)               # 仅生产帧（0-based）
frame_labels=np.array(prod_idx)         # 原始帧号，写入 CSV 便于追溯
t_ns=t_all[SKIP:]
print("total frames:",tr_all.n_frames,"| production frames:",tr.n_frames,"atoms:",tr.n_atoms)

# ---------- 1) RMSD ----------
# 以首个生产帧为参考，在蛋白骨架上做最优叠合
ref=tr.slice(0)
bb_rmsd=md.rmsd(tr,ref,atom_indices=backbone,parallel=False)
tr.superpose(ref,atom_indices=backbone)                 # 持久叠合，后续都在该坐标系
# 配体重原子 RMSD（不再二次叠合，反映蛋白对齐后配体的实际漂移）
def rmsd_noalign(coords,refc):
    d=coords-refc
    return np.sqrt((d*d).sum(axis=(1,2))/coords.shape[1])
lig_rmsd=rmsd_noalign(tr.xyz[:,lig_heavy,:],tr.xyz[0,lig_heavy,:])
rmsd_df=pd.DataFrame({"frame":frame_labels,"time_ns":t_ns,
                      "protein_backbone_RMSD_nm":bb_rmsd,"ligand_heavy_RMSD_nm":lig_rmsd})
rmsd_df.to_csv(os.path.join(ROOT,"MD_RMSD.csv"),index=False,encoding="utf-8-sig")

# ---------- 2) RMSF（蛋白，按残基） ----------
xyz=tr.xyz
mean_xyz=xyz[:,ca,:].mean(axis=0)
ca_rmsf=np.sqrt(((xyz[:,ca,:]-mean_xyz)**2).sum(axis=2)).mean(axis=0)
# 骨架平均 RMSF
bb_set=set(int(i) for i in backbone)
res_ids=[]; res_names=[]; bb_rmsf=[]
for k,a in enumerate(ca):
    res=top.atom(int(a)).residue
    bb_in_res=[i for i in backbone if top.atom(int(i)).residue==res]
    d=xyz[:,bb_in_res,:]-xyz[:,bb_in_res,:].mean(axis=0)
    bb_rmsf.append(np.sqrt((d*d).sum(axis=2)).mean())
    res_ids.append(res.resSeq); res_names.append(res.name)
rmsf_df=pd.DataFrame({"resSeq":res_ids,"residue":res_names,
                      "CA_RMSF_nm":ca_rmsf,"backbone_RMSF_nm":bb_rmsf})
rmsf_df.to_csv(os.path.join(ROOT,"MD_RMSF.csv"),index=False,encoding="utf-8-sig")

# ---------- 3) 配体-蛋白氢键 ----------
solute=tr.atom_slice(prot_atoms+lig_atoms)
lig_set=set(lig_atoms)
per_frame_count=np.zeros(tr.n_frames,dtype=int)
import collections
occ_counter=collections.Counter()
for fi in range(tr.n_frames):
    hb=md.baker_hubbard(solute.slice(fi))
    keys_this=set()
    for donor,hh,acc in hb:
        donor,hh,acc=int(donor),int(hh),int(acc)
        d_in=donor in lig_set; a_in=acc in lig_set
        if d_in == a_in:
            # 同属配体（分子内氢键）或同属蛋白（蛋白内部），不计入
            continue
        per_frame_count[fi]+=1
        prot_idx=acc if d_in else donor
        r=top.atom(prot_idx).residue
        keys_this.add((r.resSeq,r.name,"ligDonor" if d_in else "ligAcceptor"))
    for k in keys_this:
        occ_counter[k]+=1   # 同一帧同一残基只计一次
hb_df=pd.DataFrame({"frame":frame_labels,"time_ns":t_ns,
                    "n_ligand_protein_HB":per_frame_count})
hb_df.to_csv(os.path.join(ROOT,"MD_hydrogen_bonds.csv"),index=False,encoding="utf-8-sig")
occ_df=pd.DataFrame([(rn,nm,kind,c,round(c/tr.n_frames,3))
                     for (rn,nm,kind),c in occ_counter.items()],
                    columns=["resSeq","residue","ligand_role","n_frames","occupancy"])\
                     .sort_values("occupancy",ascending=False)
occ_df.to_csv(os.path.join(ROOT,"MD_hydrogen_bond_occupancy.csv"),index=False,encoding="utf-8-sig")

# ---------- 4) 配体-蛋白最小距离 ----------
prot_heavy=[i for i in prot_atoms if is_heavy(i)]
min_dist=np.zeros(tr.n_frames); com_dist=np.zeros(tr.n_frames)
lc=solute.xyz[:,lig_heavy,:].mean(axis=1)
pc=solute.xyz[:,prot_heavy,:].mean(axis=1)
for fi in range(tr.n_frames):
    a=solute.xyz[fi,lig_heavy,:]; b=solute.xyz[fi,prot_heavy,:]
    d=np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
    min_dist[fi]=d.min()
    com_dist[fi]=np.linalg.norm(lc[fi]-pc[fi])
dist_df=pd.DataFrame({"frame":frame_labels,"time_ns":t_ns,
                      "min_lig_prot_distance_nm":min_dist,"lig_prot_COM_distance_nm":com_dist})
dist_df.to_csv(os.path.join(ROOT,"MD_lig_prot_distance.csv"),index=False,encoding="utf-8-sig")

# ---------- 5) Fig8 拼版 ----------
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"]=["Microsoft YaHei","SimHei","Arial"]
fig,axs=plt.subplots(2,2,figsize=(14,9))
ax=axs[0,0]
ax.plot(t_ns,bb_rmsd,label="Protein backbone",color="#4C78A8")
ax.plot(t_ns,lig_rmsd,label="Luteolin (heavy)",color="#E15759")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (nm)"); ax.set_title("A. RMSD",loc="left",fontweight="bold")
ax.legend(frameon=False)
ax=axs[0,1]
ax.plot(res_ids,ca_rmsf,color="#59A14F")
ax.set_xlabel("Residue"); ax.set_ylabel("Cα RMSF (nm)"); ax.set_title("B. Cα RMSF",loc="left",fontweight="bold")
ax=axs[1,0]
ax.plot(t_ns,per_frame_count,color="#B07AA1")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("H-bond count")
ax.set_title("C. Ligand–protein H-bonds",loc="left",fontweight="bold")
ax=axs[1,1]
ax.plot(t_ns,min_dist,label="Min distance",color="#F28E2B")
ax.plot(t_ns,com_dist,label="COM distance",color="#76B7B2")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Distance (nm)"); ax.set_title("D. Ligand–protein distance",loc="left",fontweight="bold")
ax.legend(frameon=False)
fig.suptitle("5KIR–luteolin MD stability",fontsize=14)
plt.tight_layout(rect=[0,0,1,0.97])
figpath=os.path.join(OUTDIR,"Fig8_MD稳定性.png")
plt.savefig(figpath,dpi=300,bbox_inches="tight",facecolor="white")
print("saved",figpath)
print("mean backbone RMSD=%.3f nm, mean ligand RMSD=%.3f nm, mean HB=%.2f"%(
    bb_rmsd.mean(),lig_rmsd.mean(),per_frame_count.mean()))
