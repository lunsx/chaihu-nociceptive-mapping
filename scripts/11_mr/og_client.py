"""OpenGWAS REST client (POST-first, array params repeated in the query string).

The access token is read from the OPENGWAS_JWT environment variable, or from a
local file given by OPENGWAS_JWT_FILE. No token is committed to the repository.
"""
import os
import time

import requests

BASE = "https://api.opengwas.io/api"


def _load_token():
    tok = os.environ.get("OPENGWAS_JWT")
    if tok:
        return tok.strip()
    path = os.environ.get("OPENGWAS_JWT_FILE")
    if path and os.path.exists(path):
        return open(path).read().strip()
    raise RuntimeError(
        "OpenGWAS token missing. Set OPENGWAS_JWT or OPENGWAS_JWT_FILE.")


class OpenGWAS:
    def __init__(self):
        token = _load_token()
        self.session = requests.Session()
        self.session.trust_env = False
        self.session.headers.update({"Authorization": f"Bearer {token}"})

    def _request(self, method, endpoint, params=None, tries=6, timeout=90):
        # requests repeats list/tuple values as multi query params when doseq=True
        url = f"{BASE}/{endpoint}"
        if isinstance(params, dict):
            raw = list(params.items())
        else:
            raw = list(params or [])
        payload = [(k, v) for k, v in raw if v is not None]
        for attempt in range(tries):
            try:
                r = self.session.request(method, url, params=payload, timeout=timeout)
                if r.status_code == 429:
                    time.sleep(float(r.headers.get("retry-after", 8)))
                    continue
                if r.status_code >= 500:
                    time.sleep(2 + attempt * 2)
                    continue
                r.raise_for_status()
                return r.json()
            except requests.RequestException:
                if attempt == tries - 1:
                    raise
                time.sleep(2 + attempt * 2)

    def user(self):
        return self._request("GET", "user")

    def gwasinfo(self, ids):
        if isinstance(ids, str):
            ids = [ids]
        return self._request("POST", "gwasinfo", [("id", list(ids))])

    def tophits(self, gwas_id, pval=5e-8, clump=1, preclumped=1,
                r2=0.001, kb=5000, pop="EUR"):
        params = [
            ("id", gwas_id), ("pval", pval), ("clump", clump),
            ("preclumped", preclumped), ("r2", r2), ("kb", kb), ("pop", pop),
        ]
        return self._request("POST", "tophits", params)

    def associations(self, gwas_id, variants, proxies=0, population="EUR",
                     r2=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3):
        params = [("id", gwas_id), ("variant", list(variants)),
                  ("proxies", proxies), ("population", population),
                  ("r2", r2), ("align_alleles", align_alleles),
                  ("palindromes", palindromes), ("maf_threshold", maf_threshold)]
        return self._request("POST", "associations", params)

    def clump(self, rsids, pvals, pthresh=5e-8, r2=0.001, kb=5000, pop="EUR"):
        params = [("rsid", list(rsids)), ("pval", list(pvals)),
                  ("pthresh", pthresh), ("r2", r2), ("kb", kb), ("pop", pop)]
        return self._request("POST", "ld/clump", params)
