"""Compute SHA256 checksums for the bundled small data files.

Run:  python tests/make_checksums.py
Writes data/processed_small/CHECKSUMS.sha256 and verifies a second pass.
"""
import hashlib
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SMALL = ROOT / "data" / "processed_small"
OUT = SMALL / "CHECKSUMS.sha256"


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    lines = []
    for path in sorted(SMALL.iterdir()):
        if path.is_file() and path.name != "CHECKSUMS.sha256":
            lines.append(f"{sha256(path)}  {path.name}")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(lines)} checksums to {OUT.name}")


if __name__ == "__main__":
    main()
