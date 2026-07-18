"""Bundled reference data: local clone, user cache dir, or download.

- ``Cell_marker_Seq.xlsx`` is fetched from the CellMarker project HTTP URL.
- Other files are expected on a GitHub ``raw`` mirror (``data/`` folder).

Override locations:
- ``CELLTYPER_DATA_DIR``: use this folder (optional downloads if files missing and not offline).
- ``CELLTYPER_DATA_BASE_URL``: base URL for GitHub raw files, e.g.
  ``https://raw.githubusercontent.com/ORG/REPO/main/data`` (no trailing slash required).
- ``CELLTYPER_OFFLINE=1``: do not download; require existing files.
"""
from __future__ import annotations

import os
from pathlib import Path
from urllib.parse import urljoin

import platformdirs
import requests

CELLMARKER_SEQ_XLSX_URL = (
    "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Seq.xlsx"
)

# Replace with your released repo before publishing to PyPI, or set CELLTYPER_DATA_BASE_URL.
DEFAULT_GITHUB_DATA_BASE_URL = (
    "https://raw.githubusercontent.com/PLACEHOLDER_ORG/PLACEHOLDER_REPO/main/data"
)

# All files required for ``run_celltyper`` + CellTyper / geneRetriever flows.
REQUIRED_DATA_FILES: tuple[str, ...] = (
    "Cell_marker_Seq.xlsx",
    "PanglaoDB_markers_27_Mar_2020.tsv",
    "singleCellBase_dataset.txt",
    "proteinatlas.tsv",
    "interaction_consensus.tsv",
    "cl-basic.obo",
    "arabidopsis_thaliana.marker_fd.xlsx",
    "human_celltype.json",
    "musmusculus_celltype.json",
    "arabidopsis_celltype.json",
    "human_geneDB.csv",
    "mouse_geneDB.csv",
    "arabidopsis_geneDB.csv",
)


def data_file(filename: str) -> str:
    """Path string for a file under the active data directory (env or ``data``)."""
    base = os.environ.get("CELLTYPER_DATA_DIR", "data")
    return str(Path(base) / filename)


def _have_all(root: Path) -> bool:
    return all((root / name).is_file() and (root / name).stat().st_size > 0 for name in REQUIRED_DATA_FILES)


def _default_cache_root() -> Path:
    return Path(platformdirs.user_data_dir("celltyper")) / "data"


def _github_base_url() -> str:
    base = (os.environ.get("CELLTYPER_DATA_BASE_URL") or DEFAULT_GITHUB_DATA_BASE_URL).rstrip("/")
    if "PLACEHOLDER_ORG" in base or "PLACEHOLDER_REPO" in base:
        raise RuntimeError(
            "Set environment variable CELLTYPER_DATA_BASE_URL to your GitHub raw URL "
            "(e.g. https://raw.githubusercontent.com/ORG/REPO/main/data), or edit "
            "DEFAULT_GITHUB_DATA_BASE_URL in celltyper_data.py before release."
        )
    return base


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
    tmp.replace(dest)


def _download_missing(root: Path) -> None:
    if os.environ.get("CELLTYPER_OFFLINE", "").strip() in ("1", "true", "yes"):
        missing = [n for n in REQUIRED_DATA_FILES if not (root / n).is_file()]
        if missing:
            raise FileNotFoundError(
                "CELLTYPER_OFFLINE is set but data files are missing: "
                + ", ".join(missing)
            )
        return

    gh = _github_base_url()

    for name in REQUIRED_DATA_FILES:
        dest = root / name
        if dest.is_file() and dest.stat().st_size > 0:
            continue
        if name == "Cell_marker_Seq.xlsx":
            url = CELLMARKER_SEQ_XLSX_URL
        else:
            url = urljoin(gh + "/", name.replace("\\", "/"))
        _download(url, dest)


def ensure_celltyper_data() -> Path:
    """Ensure REQUIRED_DATA_FILES exist; set ``CELLTYPER_DATA_DIR`` and return its path."""
    explicit = os.environ.get("CELLTYPER_DATA_DIR")
    if explicit:
        root = Path(explicit).expanduser().resolve()
        root.mkdir(parents=True, exist_ok=True)
        if not _have_all(root):
            _download_missing(root)
        if not _have_all(root):
            raise FileNotFoundError(
                f"Incomplete data under CELLTYPER_DATA_DIR={root!s}. "
                "Check network, CELLTYPER_DATA_BASE_URL, and required filenames."
            )
        os.environ["CELLTYPER_DATA_DIR"] = str(root)
        return root

    cwd_data = Path("data").resolve()
    if cwd_data.is_dir() and _have_all(cwd_data):
        os.environ["CELLTYPER_DATA_DIR"] = str(cwd_data)
        return cwd_data

    root = _default_cache_root()
    root.mkdir(parents=True, exist_ok=True)
    _download_missing(root)
    if not _have_all(root):
        raise FileNotFoundError(
            f"Failed to obtain all reference files under {root}. "
            "Set CELLTYPER_DATA_BASE_URL or CELLTYPER_DATA_DIR."
        )
    os.environ["CELLTYPER_DATA_DIR"] = str(root)
    return root
