# -*- coding: utf-8 -*-
import os, json, requests, numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["SimHei"]; plt.rcParams["axes.unicode_minus"] = False

OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
PP = os.path.join(OUT, "ppi"); os.makedirs(PP, exist_ok=True)
PX = {"http": None, "https": None}

# ---- 237 formula human genes ----
M = pd.read_csv(os.path.join(OUT,"TCMSP_全方成分主表.csv"))
T = pd.read_csv(os.path.join(OUT,"TCMSP_成分靶点主表.csv"))
mp = pd.read_csv(os.path.join(OUT,"target_gene_mapping.csv")).dropna(subset=["gene"])
T["gene"] = T["target_name"].map(dict(zip(mp.target_name, mp.gene)))
fmol = set(M[M.pass_filter].MOL_ID)
genes = sorted(set(T[T.MOL_ID.isin(fmol)].gene.dropna()))
print("input genes:", len(genes))

# ---- STRING network (POST) ----
r = requests.post("https://string-db.org/api/json/network",
                  data={"identifiers": "%0d".join(genes), "species": 9606,
                        "required_score": 400, "caller_identity": "chaihu_analysis"},
                  proxies=PX, timeout=90)
r.raise_for_status(); edges = r.json()
print("STRING interactions:", len(edges))
assert len(edges) > 50, "too few interactions, abort"

G = nx.Graph()
for g in genes: G.add_node(g)
for e in edges:
    a, b = e["preferredName_A"], e["preferredName_B"]
    G.add_edge(a, b, score=float(e["score"]))
print("graph nodes/edges:", G.number_of_nodes(), G.number_of_edges())

# ---- topology ----
deg = dict(G.degree())
wdeg = dict(G.degree(weight="score"))
bet = nx.betweenness_centrality(G)
close = nx.closeness_centrality(G)
eig = nx.eigenvector_centrality(G, max_iter=1000)
nodes = pd.DataFrame({
    "gene": list(G.nodes()),
    "degree": [deg[n] for n in G.nodes],
    "weighted_degree": [round(wdeg[n],3) for n in G.nodes],
    "betweenness": [round(bet[n],4) for n in G.nodes],
    "closeness": [round(close[n],4) for n in G.nodes],
    "eigenvector": [round(eig[n],4) for n in G.nodes],
})
# composite rank (mean of percentile ranks)
for c in ["degree","betweenness","closeness","eigenvector"]:
    nodes[c+"_pct"] = nodes[c].rank(pct=True)
nodes["hub_score"] = nodes[[c+"_pct" for c in ["degree","betweenness","closeness","eigenvector"]]].mean(axis=1)
nodes = nodes.sort_values("hub_score", ascending=False)
# mark pain set
PAIN = set(["TRPV1","PTGS2","TNF","IL1B","IL6","IL1A","PTGES","PTGER3","PTGS1","TRPA1",
            "TACR1","CALCRL","NGF","SCN9A","SCN5A","SCN10A","TRPV3","TRPC3","IL1R1"])
nodes["is_pain_gene"] = nodes.gene.isin(PAIN)
nodes.to_csv(os.path.join(PP,"PPI_节点拓扑.csv"), index=False, encoding="utf-8-sig")

edge_df = pd.DataFrame([{"A":a,"B":b,"score":d["score"]} for a,b,d in G.edges(data=True)])
edge_df.to_csv(os.path.join(PP,"PPI_边表.csv"), index=False, encoding="utf-8-sig")

# ---- plot: top 60 by degree for readability ----
topN = nodes.head(60)["gene"].tolist()
sub = G.subgraph(topN)
pos = nx.spring_layout(sub, seed=7, k=0.55)
fig, ax = plt.subplots(figsize=(14,11)); ax.axis("off")
dmax = max(deg[n] for n in topN)
sizes = [300 + 1700*deg[n]/dmax for n in topN]
cols = ["#C0392B" if n in PAIN else "#4A90D9" for n in topN]
nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.25, edge_color="#7F8C8D")
nx.draw_networkx_nodes(sub, pos, ax=ax, nodelist=topN, node_size=sizes,
                       node_color=cols, alpha=0.9)
# white bold labels on red nodes, normal on blue
pain_in = [n for n in topN if n in PAIN]; other_in = [n for n in topN if n not in PAIN]
nx.draw_networkx_labels(sub, pos, labels={n:n for n in pain_in}, ax=ax,
                        font_size=10, font_weight="bold", font_color="white")
nx.draw_networkx_labels(sub, pos, labels={n:n for n in other_in}, ax=ax, font_size=9)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0],marker="o",color="w",markerfacecolor="#C0392B",markersize=13,label="痛觉/炎症基因"),
                   Line2D([0],[0],marker="o",color="w",markerfacecolor="#4A90D9",markersize=13,label="其他成分靶点")],
          loc="lower left", fontsize=11, frameon=True)
ax.set_title("成分靶点 PPI 网络（Top 60 核心节点；红色=痛觉/炎症基因）", fontsize=14)
plt.tight_layout(); plt.savefig(os.path.join(PP,"fig_PPI网络.png"), dpi=200); plt.close()

print("\nTop 20 hubs:")
print(nodes.head(20)[["gene","degree","betweenness","eigenvector","hub_score","is_pain_gene"]].to_string(index=False))
print("\nfiles:", os.listdir(PP))
