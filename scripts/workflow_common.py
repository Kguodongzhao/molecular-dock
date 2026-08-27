from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = PACKAGE_ROOT / "config" / "docking.local.json"


def load_config(path: str | Path | None = None) -> dict:
    config_path = Path(path).resolve() if path else DEFAULT_CONFIG

    if not config_path.is_file():
        raise FileNotFoundError(
            f"Config file not found: {config_path}\n"
            "Copy config/docking.example.json to "
            "config/docking.local.json and edit it, "
            "or pass --config PATH."
        )

    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["_config_path"] = str(config_path)
    return config


def package_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PACKAGE_ROOT / path


def find_vina(config: dict) -> Path:
    override = os.environ.get("VINA_EXE")

    if override:
        override_path = Path(override).expanduser()
        if override_path.is_file():
            return override_path

        raise FileNotFoundError(
            f"VINA_EXE does not point to a file: {override_path}"
        )

    configured_path = config.get("vina_relative_path")
    if configured_path:
        bundled = package_path(configured_path)
        if bundled.is_file():
            return bundled

    executable = shutil.which("vina") or shutil.which("vina.exe")
    if executable:
        return Path(executable)

    raise FileNotFoundError(
        "Cannot find AutoDock Vina. Set VINA_EXE, add vina to PATH, "
        "or set vina_relative_path in the configuration file."
    )


def find_command(name: str) -> str:
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(
            f"Cannot find command '{name}'. Activate the molecular-docking conda environment."
        )
    return found


def run_command(args: list[str], timeout: int = 3600) -> subprocess.CompletedProcess:
    return subprocess.run(
        args,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


def numeric_key(path: Path) -> tuple[int, int | str]:
    try:
        return 0, int(path.stem.split("_")[0])
    except ValueError:
        return 1, path.stem


def parse_case_list(value: str | None) -> set[str] | None:
    if not value:
        return None
    return {item.strip().replace("\\", "/") for item in value.split(",") if item.strip()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdb_atom_lines(path: Path) -> list[str]:
    return [
        line
        for line in path.read_text(encoding="ascii", errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54
    ]


def pdb_elements(path: Path) -> set[str]:
    elements = set()
    for line in pdb_atom_lines(path):
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            atom_name = re.sub(r"[^A-Za-z]", "", line[12:16]).upper()
            element = atom_name[:2] if atom_name[:2] in {"BE", "CL", "BR"} else atom_name[:1]
        if element:
            elements.add(element)
    return elements


def pdbqt_torsdof(path: Path) -> int | None:
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("TORSDOF"):
            try:
                return int(line.split()[1])
            except (IndexError, ValueError):
                return None
    return None


def vina_scores(path: Path) -> list[float]:
    scores = []
    for line in path.read_text(encoding="ascii", errors="replace").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            try:
                scores.append(float(line.split()[3]))
            except (IndexError, ValueError):
                return []
    return scores


def pdbqt_atom_count(path: Path) -> int:
    return sum(
        line.startswith(("ATOM  ", "HETATM"))
        for line in path.read_text(encoding="ascii", errors="replace").splitlines()
    )


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def build_complex(display_receptor: Path, ligand_pose: Path, output: Path) -> None:
    receptor = pdb_atom_lines(display_receptor)
    ligand = []
    for line in pdb_atom_lines(ligand_pose):
        ligand.append("HETATM" + line[6:] if line.startswith("ATOM  ") else line)
    output.write_text(
        "\n".join(receptor + ["TER"] + ligand + ["END", ""]),
        encoding="ascii",
    )
