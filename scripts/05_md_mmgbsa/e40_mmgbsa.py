# -*- coding: utf-8 -*-
"""
e40  单轨迹 MM-GBSA（5KIR-luteolin），OBC2，NoCutoff，省略熵
公式：dG_bind = G_complex - G_protein - G_ligand（同一帧坐标，不做分别最小化）
产出：MMGBSA_per_frame.csv / MMGBSA_summary.csv
依赖：OpenMM；读取 system.xml、modeller.pkl、traj.dcd
说明：
  - 蛋白力场 ff14SB、配体 OpenFF Sage 沿用原系统参数；
  - GB 半径按元素取 OBC2/mbondi 常用值（C 0.17, O 0.15, N 0.155, S 0.20,
    H 连 O/N 0.10 否则 0.12，单位 nm），属近似，已在局限中声明；
  - 非极性项用 GBSAOBCForce 的 SurfaceAreaEnergy；
  - 命令行：python e40_mmgbsa.py [stride]，默认 2（每隔一帧取一帧）。
"""
import os, sys, pickle
import numpy as np, pandas as pd
import mdtraj as md
from openmm import XmlSerializer, System, Platform, Vec3, Context, VerletIntegrator
from openmm import HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce, NonbondedForce, GBSAOBCForce
from openmm import unit
CH=unit.elementary_charge; CH2=unit.elementary_charge**2; NM=unit.nanometer; KJ=unit.kilojoule_per_mole

MD=r"C:\dockwork\md"
ROOT=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"

with open(os.path.join(MD,"system.xml")) as f:
    FULL=XmlSerializer.deserialize(f.read())
_m=pickle.load(open(os.path.join(MD,"modeller.pkl"),"rb"))
TOP=_m["topology"]; N_PROT=int(_m["n_prot"]); LIG=list(_m["lig_idx"])
ELEMS=[a.element.symbol for a in list(TOP.atoms())]

_bondmap={}
for _fc in FULL.getForces():
    if isinstance(_fc,HarmonicBondForce):
        for _i in range(_fc.getNumBonds()):
            _a,_b,_,_=_fc.getBondParameters(_i); _a=int(_a); _b=int(_b)
            _bondmap.setdefault(_a,[]).append(_b); _bondmap.setdefault(_b,[]).append(_a)

def radius(i):
    e=ELEMS[i]
    if e=="C": return 0.17
    if e=="O": return 0.15
    if e=="N": return 0.155
    if e=="S": return 0.20
    if e=="H":
        heavy=[j for j in _bondmap.get(i,[]) if ELEMS[j]!="H"]
        if heavy and ELEMS[heavy[0]] in ("O","N"): return 0.10
        return 0.12
    return 0.15

def build_subsystem(keep):
    keep=set(keep); old2new={old:k for k,old in enumerate(sorted(keep))}
    nb_src=next(f for f in FULL.getForces() if isinstance(f,NonbondedForce))
    charges={i:nb_src.getParticleParameters(i)[0].value_in_unit(CH)
             for i in range(nb_src.getNumParticles())}
    sysn=System()
    for old in sorted(keep):
        sysn.addParticle(FULL.getParticleMass(old))
    bf=HarmonicBondForce(); bf.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc,HarmonicBondForce):
            for i in range(fc.getNumBonds()):
                a,b,L,k=fc.getBondParameters(i); a=int(a); b=int(b)
                if a in keep and b in keep: bf.addBond(old2new[a],old2new[b],L,k)
    af=HarmonicAngleForce(); af.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc,HarmonicAngleForce):
            for i in range(fc.getNumAngles()):
                a,b,c,ang,k=fc.getAngleParameters(i); a,b,c=int(a),int(b),int(c)
                if a in keep and b in keep and c in keep:
                    af.addAngle(old2new[a],old2new[b],old2new[c],ang,k)
    tf=PeriodicTorsionForce(); tf.setForceGroup(1)
    for fc in FULL.getForces():
        if isinstance(fc,PeriodicTorsionForce):
            for i in range(fc.getNumTorsions()):
                a,b,c,d,per,ph,k=fc.getTorsionParameters(i); a,b,c,d=int(a),int(b),int(c),int(d)
                if a in keep and b in keep and c in keep and d in keep:
                    tf.addTorsion(old2new[a],old2new[b],old2new[c],old2new[d],per,ph,k)
    nb=NonbondedForce(); nb.setNonbondedMethod(NonbondedForce.NoCutoff); nb.setForceGroup(2)
    for old in sorted(keep):
        q,sig,eps=nb_src.getParticleParameters(old)
        nb.addParticle(q.value_in_unit(CH)*CH, sig.value_in_unit(NM)*NM, eps.value_in_unit(KJ)*KJ)
    for j in range(nb_src.getNumExceptions()):
        a,b,q,sig,eps=nb_src.getExceptionParameters(j); a=int(a); b=int(b)
        if a in keep and b in keep:
            nb.addException(old2new[a],old2new[b],
                            q.value_in_unit(CH2)*CH2, sig.value_in_unit(NM)*NM, eps.value_in_unit(KJ)*KJ)
    gb=GBSAOBCForce(); gb.setNonbondedMethod(GBSAOBCForce.NoCutoff)
    gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0); gb.setSurfaceAreaEnergy(True)
    gb.setForceGroup(3)
    for old in sorted(keep):
        gb.addParticle(charges[old],radius(old),0.7)
    for f in (bf,af,tf,nb,gb): sysn.addForce(f)
    return sysn, sorted(keep)

# DCD 前 6 帧（step 10000..60000）是 e38 平衡/约束阶段写入，生产帧从索引 6（step 70000）起
SKIP_FRAMES=6

def main(stride=2):
    species={"complex":list(range(0,LIG[-1]+1)),
             "protein":list(range(0,N_PROT)),
             "ligand":LIG}
    subs={}
    for name,keep in species.items():
        sysn,kept=build_subsystem(keep)
        ctx=Context(sysn,VerletIntegrator(0.001),Platform.getPlatformByName("CPU"))
        subs[name]=(ctx,kept); print(name,"atoms:",len(kept))
    mtop=md.Topology.from_openmm(TOP)
    tr=md.load(os.path.join(MD,"traj.dcd"),top=mtop)
    frames=list(range(SKIP_FRAMES,tr.n_frames,stride))
    print("total frames:",tr.n_frames,"| production frames evaluated:",len(frames))
    rows=[]
    for fi in frames:
        rec={"frame":fi,"time_ns":round(fi*0.04,3)}
        for name,(ctx,kept) in subs.items():
            ctx.setPositions([Vec3(*tr.xyz[fi,a]) for a in kept])
            e_bond=ctx.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(KJ)
            e_nb=ctx.getState(getEnergy=True,groups={2}).getPotentialEnergy().value_in_unit(KJ)
            e_gb=ctx.getState(getEnergy=True,groups={3}).getPotentialEnergy().value_in_unit(KJ)
            rec[name+"_bonded"]=e_bond; rec[name+"_MM"]=e_nb; rec[name+"_GBSA"]=e_gb
            rec[name+"_total"]=e_bond+e_nb+e_gb
        for comp in ["bonded","MM","GBSA","total"]:
            d=rec["complex_"+comp]-rec["protein_"+comp]-rec["ligand_"+comp]
            rec["dG_"+comp+"_kJ"]=d; rec["dG_"+comp+"_kcal"]=d/4.184
        rows.append(rec)
    pf=pd.DataFrame(rows)
    pf.to_csv(os.path.join(ROOT,"MMGBSA_per_frame.csv"),index=False,encoding="utf-8-sig")
    summ=[]
    for comp in ["bonded","MM","GBSA","total"]:
        x=pf["dG_"+comp+"_kcal"]
        summ.append({"component":comp,"mean_kcal_mol":x.mean(),"SD":x.std(ddof=1),
                     "SEM":x.std(ddof=1)/np.sqrt(len(x)),"median":x.median()})
    sf=pd.DataFrame(summ)
    sf.to_csv(os.path.join(ROOT,"MMGBSA_summary.csv"),index=False,encoding="utf-8-sig")
    print(sf.to_string(index=False)); print("n evaluated frames:",len(frames))

if __name__=="__main__":
    main(int(sys.argv[1]) if len(sys.argv)>1 else 2)
