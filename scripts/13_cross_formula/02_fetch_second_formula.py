# -*- coding: utf-8 -*-
"""Fetch the second-formula (Jinlingzi San) herb pages from TCMSP.

The two herbs share no herbs with Chaihu Shugan San. This script uses the
English registered names directly with the TCMSP detail endpoint, so it needs
no Chinese input. Raw HTML is written under the project work root.

Run from the repository; the shared config resolves the project root and works
on both the Windows workstation and the Linux host.
"""
import re
import sys
import time
from pathlib import Path

import requests

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
from config import PATHS  # noqa: E402

BASE_URL = "https://old.tcmsp-e.com/"
HERBS = {
    "Toosendan Fructus": "Toosendan_Fructus.html",
    "Corydalis Rhizoma": "Corydalis_Rhizoma.html",
}
RAW_DIR = PATHS["project_root"] / "xiaoyao_case" / "jinling_raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.proxies = {"http": None, "https": None}
session.headers.update({"User-Agent": "Mozilla/5.0"})
token = re.search(r'name="token" value="([0-9a-f]+)"',
                  session.get(BASE_URL + "tcmsp.php", timeout=25).text).group(1)


def fetch_detail(english_name):
    for _ in range(3):
        resp = session.get(
            BASE_URL + "tcmspsearch.php",
            params={"qr": english_name, "qsr": "herb_en_name", "token": token},
            timeout=90,
        )
        if resp.status_code == 200 and "platform_table" not in resp.text and len(resp.text) > 50000:
            return resp.text
        time.sleep(2)
    return None


for name, out_name in HERBS.items():
    html = fetch_detail(name)
    if html is None:
        print(name, "FAILED")
        continue
    (RAW_DIR / out_name).write_text(html, encoding="utf-8")
    print(name, "saved", len(html) // 1024, "KB")
