from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


EXPECTED_CONDA_ENV = "molecular-docking"
PREPROCESSING_ROOT = Path(__file__).resolve().parent


def require_molecular_docking_environment() -> None:
    active = os.environ.get("CONDA_DEFAULT_ENV", "")
    if active == EXPECTED_CONDA_ENV or Path(active).name == EXPECTED_CONDA_ENV:
        return
    raise RuntimeError(
        "This script must run in the 'molecular-docking' conda environment. "
        "Use preprocessing/run.ps1 or activate the environment first."
    )


def find_command(name: str) -> str:
    command = shutil.which(name)
    if not command:
        raise FileNotFoundError(
            f"Cannot find '{name}' in the molecular-docking environment."
        )
    return command


def run_command(
    command: list[str], *, timeout: int = 3600
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


def command_error(completed: subprocess.CompletedProcess[str]) -> str:
    return (completed.stderr or completed.stdout).strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_json(path: Path, value: dict) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def safe_stem(value: str, fallback: str = "molecule") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    return cleaned or fallback


def pdbqt_atom_count(path: Path) -> int:
    return sum(
        line.startswith(("ATOM  ", "HETATM"))
        for line in path.read_text(encoding="ascii", errors="replace").splitlines()
    )


def pdbqt_torsdof(path: Path) -> int | None:
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("TORSDOF"):
            try:
                return int(line.split()[1])
            except (IndexError, ValueError):
                return None
    return None


def tool_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    try:
        from importlib.metadata import version

        versions["meeko"] = version("meeko")
        versions["rdkit"] = version("rdkit")
    except Exception:
        pass
    obabel = shutil.which("obabel")
    if obabel:
        completed = run_command([obabel, "-V"], timeout=30)
        versions["openbabel"] = (completed.stdout or completed.stderr).strip()
    return versions
