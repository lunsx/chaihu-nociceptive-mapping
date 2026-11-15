"""Shared path and environment configuration.

Import at the top of an analysis script:

    from config import PATHS

Environment variables override the defaults, so the same script can run on the
Windows workstation, the Linux cloud host, or a fresh clone without editing:

    CHAIHU_PROJECT_ROOT   project work root (contains the analysis folders)
    CHAIHU_DOCK_ROOT      docking / heavy-IO root (default C:\\dockwork)
    CHAIHU_VINA           AutoDock Vina executable
"""
import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _default_project_root() -> Path:
    # repository/scripts/config.py: parents[1] is the repo, parents[2] is the
    # project work root that holds the analysis folders and this repository.
    here = Path(__file__).resolve()
    candidates = [
        here.parents[2],
        Path(r"C:\Users\Administrator\Desktop\草药\新选题_醋香附胁痛\计算分析"),
        Path("/home/user/compute_work"),
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    return here.parents[2]


PROJECT_ROOT = Path(os.environ.get("CHAIHU_PROJECT_ROOT", _default_project_root()))
DOCK_ROOT = Path(os.environ.get("CHAIHU_DOCK_ROOT", r"C:\dockwork"))

PATHS = {
    "repo_root": REPO_ROOT,
    "project_root": PROJECT_ROOT,
    "dock_root": DOCK_ROOT,
    "dude_bench": DOCK_ROOT / "dude_bench",
    "masld": DOCK_ROOT / "masld",
    "psc": DOCK_ROOT / "psc",
    "md": DOCK_ROOT / "md",
    "healthy_liver": PROJECT_ROOT / "healthy_liver",
    "singlecell": PROJECT_ROOT / "singlecell",
    "vina": Path(os.environ.get("CHAIHU_VINA", DOCK_ROOT / "vina.exe")),
    "data": REPO_ROOT / "data",
    "processed_small": REPO_ROOT / "data" / "processed_small",
    "figures": REPO_ROOT / "figures",
}


def ensure_dir(key: str) -> Path:
    path = PATHS[key]
    path.mkdir(parents=True, exist_ok=True)
    return path


if __name__ == "__main__":
    for key, value in PATHS.items():
        print(f"{key:16} {value}  exists={value.exists()}")
