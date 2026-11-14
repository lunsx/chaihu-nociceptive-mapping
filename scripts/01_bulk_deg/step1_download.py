# -*- coding: utf-8 -*-
import os, urllib.request, ssl, gzip
OUT = r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
URLS = {
 "GSE89632_series_matrix.txt.gz":
   "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89632/matrix/GSE89632_series_matrix.txt.gz",
 "GPL14951_family.soft.gz":
   "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL14nnn/GPL14951/annot/GPL14951.annot.gz",
}
def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
    op = urllib.request.build_opener(urllib.request.ProxyHandler({}),
                                     urllib.request.HTTPSHandler(context=ctx))
    with op.open(req, timeout=120) as r:
        return r.read()
for name,url in URLS.items():
    dst=os.path.join(OUT,name)
    try:
        data=fetch(url)
        open(dst,"wb").write(data)
        print(f"[OK] {name} {len(data)//1024}KB")
    except Exception as e:
        print(f"[ERR] {name} {e}")
