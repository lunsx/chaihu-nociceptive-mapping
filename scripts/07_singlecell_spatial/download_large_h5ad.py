# -*- coding: utf-8 -*-
"""Download the large source h5ad files listed in large_files_manifest.csv.

The large h5ad files are not redistributed in the repository. This script reads
the manifest and downloads each file to its config location. An optional proxy
is honored through the standard HTTPS_PROXY environment variable. No paths are
hard-coded.

Run:
  python scripts/07_singlecell_spatial/download_large_h5ad
  $env:HTTPS_PROXY="http://127.0.0.1:8902"; python scripts/07_singlecell_spatial/download_large_h5ad
"""
import csv
import os
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import PATHS, REPO_ROOT  # noqa: E402

MANIFEST = REPO_ROOT / "data" / "raw_external_links" / "large_files_manifest.csv"

# Map a config location expression to a concrete directory.
KEY_TO_DIR = {
    "healthy_liver": PATHS["healthy_liver"],
    "singlecell": PATHS["singlecell"],
    "psc": PATHS["psc"],
}


def resolve_target(config_location):
    # config_location looks like PATHS["key"]/file.h5ad
    token = config_location.split("]", 1)[0].split('"')[1]
    name = config_location.split("/", 1)[1]
    return KEY_TO_DIR[token] / name


def download(url, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(target, "wb") as fh:
        total = int(r.headers.get("Content-Length", "0"))
        got = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
            got += len(chunk)
            if total:
                print(f"  {got/1e6:.0f}/{total/1e6:.0f} MB", flush=True)
    print("saved", target, target.stat().st_size, "bytes")


def main():
    with open(MANIFEST, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        target = resolve_target(row["config_location"])
        if target.exists() and target.stat().st_size == int(row["size_bytes"]):
            print("already present:", target)
            continue
        print("GET", row["source_url"])
        download(row["source_url"], target)


if __name__ == "__main__":
    main()
