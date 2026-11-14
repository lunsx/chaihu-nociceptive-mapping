# -*- coding: utf-8 -*-
import os, requests, numpy as np, json
PX={"http":None,"https":None}
DOCK=r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析\docking"
LIG=os.path.join(DOCK,"ligands"); os.makedirs(LIG,exist_ok=True)
def name2cid(name):
    r=requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON",
                   proxies=PX,timeout=40)
    return r.json()["IdentifierList"]["CID"][0]
ligs={"quercetin":5280343,"luteolin":5280445,"kaempferol":5280802,"paeoniflorin":442574}
for n in ["cyperotundone","alpha-cyperone"]:
    try: ligs[n]=name2cid(n)
    except Exception as e: print("name fail",n,e)
print("ligand CIDs:",ligs)
for n,cid in ligs.items():
    p=os.path.join(LIG,n+".sdf")
    r=requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d",
                   proxies=PX,timeout=60)
    if r.status_code==200 and len(r.content)>200:
        open(p,"wb").write(r.content); print(n,cid,"SDF",len(r.content))
    else:
        print(n,cid,"no 3D, status",r.status_code)

# pocket centers from co-crystal ligands
def center(pdb,hetid):
    xs=[];ys=[];zs=[]
    for ln in open(pdb):
        if ln.startswith("HETATM") and ln[17:20].strip()==hetid:
            xs.append(float(ln[30:38]));ys.append(float(ln[38:46]));zs.append(float(ln[46:54]))
    return (float(np.mean(xs)),float(np.mean(ys)),float(np.mean(zs))),len(xs)
c2az5,n2=center(os.path.join(DOCK,"receptors","2AZ5.pdb"),"307")
c5kir,n5=center(os.path.join(DOCK,"receptors","5KIR.pdb"),"RCX")
print("\n2AZ5 ligand307 center:",c2az5,"atoms",n2)
print("5KIR RCX center:",c5kir,"atoms",n5)
json.dump({"2AZ5":c2az5,"5KIR":c5kir},open(os.path.join(DOCK,"pocket_centers.json"),"w"),indent=1)
