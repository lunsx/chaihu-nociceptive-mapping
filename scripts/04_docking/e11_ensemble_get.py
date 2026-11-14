# -*- coding: utf-8 -*-
import os, time, requests, collections
PX = {"http": None, "https": None}
work = r"C:\dockwork\ensemble"
os.makedirs(work, exist_ok=True)
PID = "6COX"

def get(url, tries=8):
    for i in range(tries):
        try:
            r = requests.get(url, proxies=PX, timeout=60)
            if r.status_code == 200:
                return r
        except Exception:
            pass
        time.sleep(1.2 * (i + 1))
    raise RuntimeError("download fail " + url)

r = get(f"https://files.rcsb.org/download/{PID}.pdb")
pdb = os.path.join(work, PID + ".pdb")
with open(pdb, "w", encoding="utf-8", newline="\n") as f:
    f.write(r.text)
print("downloaded", PID, "chars:", len(r.text))

groups = collections.Counter()
for line in r.text.splitlines():
    if line.startswith("HETATM") and line[17:20].strip() != "HOH":
        groups[(line[17:20].strip(), line[22:26].strip(), line[21].strip())] += 1
print("non-water HET groups:")
for k, v in groups.items():
    print(" ", k, "atoms:", v)
