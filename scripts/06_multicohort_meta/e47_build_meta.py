# -*- coding: utf-8 -*-
"""解析各队列 series matrix，建立“表达列样本ID -> 分组/分级”元数据并缓存。"""
import gzip, json, os, re
import requests

PROXIES = {"http": "http://127.0.0.1:8902", "https": "http://127.0.0.1:8902"}
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"
META_DIR = os.path.join("multicohort", "_meta")
os.makedirs(META_DIR, exist_ok=True)

COHORTS = ["GSE281797", "GSE260666", "GSE147304", "GSE193066", "GSE192959"]

def matrix_url(gse):
    return f"{GEO}/{gse[:-3]}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz"

def get_matrix(gse):
    local = os.path.join(META_DIR, f"{gse}_series_matrix.txt.gz")
    if not os.path.exists(local):
        r = requests.get(matrix_url(gse), proxies=PROXIES, timeout=300)
        r.raise_for_status()
        with open(local, "wb") as f:
            f.write(r.content)
    with gzip.open(local, "rt", encoding="utf-8", errors="replace") as f:
        return f.read().read() if False else f.read()

def parse_series(txt):
    """返回每个 GSM 的 dict：title + 所有 characteristics。"""
    lines = txt.splitlines()
    gsm, titles = [], []
    char_blocks = []
    for ln in lines:
        if ln.startswith("!Sample_geo_accession"):
            gsm = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            titles = [x.strip('"') for x in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            char_blocks.append([x.strip('"') for x in ln.split("\t")[1:]])
    out = {}
    for i, g in enumerate(gsm):
        d = {"title": titles[i] if i < len(titles) else ""}
        for block in char_blocks:
            if i < len(block):
                v = block[i]
                if ":" in v:
                    k, val = v.split(":", 1)
                    d[k.strip()] = val.strip()
                else:
                    d[f"ch_{len(d)}"] = v
        out[g] = d
    return out

def num(v):
    try:
        return int(v)
    except Exception:
        try:
            return float(v)
        except Exception:
            return None

all_meta = {}

# GSE281797: title HumanLiver_S308 -> key S308
m = parse_series(get_matrix("GSE281797"))
d = {}
for g, info in m.items():
    mm = re.search(r"(S\d+)", info.get("title", ""))
    if not mm:
        continue
    d[mm.group(1)] = {
        "gsm": g,
        "diagnosis": info.get("diagnosis", ""),
        "inflammation": num(info.get("inflammation grade")),
        "fibrosis": num(info.get("fibrosis grade")),
        "steatosis": num(info.get("steatosis grade")),
        "ballooning": num(info.get("ballooning grade")),
        "nas": num(info.get("nafld activity_score")),
    }
all_meta["GSE281797"] = d

# GSE260666: title Control_rep1 [21120D1] -> key 21120D1
m = parse_series(get_matrix("GSE260666"))
d = {}
state_map = {
    "healthy control": "Control",
    "non-alcoholic fatty liver disease (NAFLD)": "NAFL",
    "non-alcoholic steatohepatitis (NASH)": "NASH",
}
for g, info in m.items():
    mm = re.search(r"\[([^\]]+)\]", info.get("title", ""))
    if not mm:
        continue
    d[mm.group(1)] = {
        "gsm": g,
        "group": state_map.get(info.get("disease state", ""), info.get("disease state", "")),
    }
all_meta["GSE260666"] = d

# GSE147304: 直接按表达列前缀分类（NASH* vs 其余），这里存 GSM 级，分析时按列名
m = parse_series(get_matrix("GSE147304"))
d = {}
for g, info in m.items():
    d[g] = {"title": info.get("title", ""), "nash": info.get("with nash", "")}
all_meta["GSE147304"] = d

# GSE193066: title HUnafld001 -> key HUnafld001
m = parse_series(get_matrix("GSE193066"))
d = {}
for g, info in m.items():
    key = info.get("title", "").strip()
    d[key] = {
        "gsm": g,
        "fibrosis": num(info.get("fibrosis stage")),
        "nas": num(info.get("nafld activity score")),
    }
all_meta["GSE193066"] = d

# GSE192959: title TUnafld001
m = parse_series(get_matrix("GSE192959"))
d = {}
for g, info in m.items():
    key = info.get("title", "").strip()
    d[key] = {
        "gsm": g,
        "fibrosis": num(info.get("fibrosis stage")),
        "nas": num(info.get("nafld activity score")),
    }
all_meta["GSE192959"] = d

with open(os.path.join(META_DIR, "cohort_metadata.json"), "w", encoding="utf-8") as f:
    json.dump(all_meta, f, ensure_ascii=False, indent=2)

for c, d in all_meta.items():
    print("=" * 80)
    print(c, "samples:", len(d))
    keys = list(d)[:5]
    for k in keys:
        print(" ", k, d[k])
