# -*- coding: utf-8 -*-
import os, requests, collections
PX={"http":None,"https":None}
DOCK=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking"
REC=os.path.join(DOCK,"receptors"); os.makedirs(REC,exist_ok=True)
pdbs={"1ITB":"IL1受体","2AZ5":"TNF","5KIR":"PTGS2","8GFA":"TRPV1"}
for pid,desc in pdbs.items():
    p=os.path.join(REC,pid+".pdb")
    if not os.path.exists(p):
        r=requests.get(f"https://files.rcsb.org/download/{pid}.pdb",proxies=PX,timeout=90)
        open(p,"wb").write(r.content)
    lines=open(p).read().splitlines()
    chains=set(); het=collections.Counter(); hetnames={}
    for ln in lines:
        if ln.startswith("ATOM"): chains.add(ln[21])
        if ln.startswith("HETATM"):
            res=ln[17:20].strip(); het[res]+=1
            hetnames.setdefault(res,ln[17:20].strip())
    print(f"{pid} ({desc}): {os.path.getsize(p)} bytes, chains={sorted(chains)}")
    print("  HETATM residues:",dict(het))
