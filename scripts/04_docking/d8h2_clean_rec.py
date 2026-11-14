# -*- coding: utf-8 -*-
import os
DOCK=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking"
REC=os.path.join(DOCK,"receptors")
plans={"2AZ5":["A","B","C"],"5KIR":["A"],"1ITB":["A","B"]}
for pid,chains in plans.items():
    src=os.path.join(REC,pid+".pdb")
    dst=os.path.join(REC,pid+"_clean.pdb")
    kept=0; nres=set()
    fh=open(src); oh=open(dst,"w")
    for ln in fh:
        if ln.startswith("ATOM") and ln[21] in chains:
            oh.write(ln); kept+=1; nres.add((ln[21],ln[22:26].strip()))
    oh.close(); fh.close()
    print(pid,"kept_atoms",kept,"residues",len(nres))
