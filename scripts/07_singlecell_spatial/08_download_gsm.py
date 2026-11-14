# -*- coding: utf-8 -*-
"""Download the remaining GSE344087 per-sample h5 files (S008-S019).

Per-sample files are only 8-14 MB, which avoids the 170 MB RAW.tar and the
partial-tar extraction used earlier. Resumable and idempotent.

Run:  python download_gsm.py
"""
import pathlib
import sys
import time

import requests

DST = pathlib.Path(r"C:\dockwork\masld\gsm")
DST.mkdir(parents=True, exist_ok=True)
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9970nnn/GSM{gsm}/suppl/{name}"
SAMPLES = [
    (9970283, "MASLD_S008"), (9970284, "MASLD_S009"), (9970285, "MASLD_S010"),
    (9970286, "MASLD_S011"), (9970287, "MASLD_S012"), (9970288, "MASLD_S013"),
    (9970289, "MASLD_S014"), (9970290, "MASLD_S015"), (9970291, "MASLD_S016"),
    (9970292, "MASLD_S017"), (9970293, "MASLD_S018"), (9970294, "MASLD_S019"),
]
S = requests.Session()
S.trust_env = True   # local proxy 127.0.0.1:8902 via env


def fetch(gsm, tag):
    name = f"GSM{gsm}_{tag}_LIVER_filtered_feature_bc_matrix.h5"
    url = BASE.format(gsm=gsm, name=name)
    out = DST / f"GSM{gsm}_{tag}.h5"
    if out.exists() and out.stat().st_size > 1_000_000:
        print(f"skip {out.name} ({out.stat().st_size} bytes)", flush=True)
        return True
    t0 = time.time()
    with S.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = out.with_suffix(".part")
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(65536):
                fh.write(chunk)
        tmp.replace(out)
    dt = time.time() - t0
    print(f"{out.name}: {out.stat().st_size} bytes in {dt:.0f}s "
          f"({out.stat().st_size/1e6/max(dt,1e-9):.2f} MB/s)", flush=True)
    return True


def main():
    ok = 0
    for gsm, tag in SAMPLES:
        try:
            if fetch(gsm, tag):
                ok += 1
        except Exception as exc:  # noqa: BLE001  network flakiness is expected
            print(f"FAIL GSM{gsm}: {exc}", file=sys.stderr, flush=True)
    print(f"done {ok}/{len(SAMPLES)}", flush=True)


if __name__ == "__main__":
    main()
