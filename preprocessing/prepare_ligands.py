from __future__ import annotations

import argparse
import concurrent.futures
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from preprocessing_common import (
    PREPROCESSING_ROOT,
    command_error,
    find_command,
    pdbqt_atom_count,
    pdbqt_torsdof,
    require_molecular_docking_environment,
    run_command,
    safe_stem,
    sha256,
    tool_versions,
    write_json,
    write_tsv,
)


SUPPORTED_EXTENSIONS = {
    ".sdf",
    ".sd",
    ".mol",
    ".mol2",
    ".pdb",
    ".smi",
    ".smiles",
    ".inchi",
}
FIELDS = [
    "source",
    "molecule_index",
    "name",
    "status",
    "output_pdbqt",
    "torsdof",
    "atom_count",
    "rigidified",
    "elapsed_seconds",
    "message",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare raw small molecules as AutoDock Vina PDBQT ligands."
    )
    parser.add_argument("--input", required=True, help="Input molecule file or directory")
    parser.add_argument(
        "--output-root",
        default=str(PREPROCESSING_ROOT / "output" / "ligands_pdbqt"),
    )
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--optimize-steps", type=int, default=500)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--no-recursive", action="store_true")
    parser.add_argument("--keep-salts", action="store_true")
    parser.add_argument("--keep-existing-3d", action="store_true")
    parser.add_argument("--rigid-macrocycles", action="store_true")
    parser.add_argument(
        "--rigid-if-over",
        type=int,
        metavar="N",
        help="Rigidify a generated ligand when TORSDOF is greater than N.",
    )
    parser.add_argument(
        "--charge-model",
        choices=("gasteiger", "zero"),
        default="gasteiger",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-intermediate", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1")
    if args.optimize_steps < 0:
        parser.error("--optimize-steps cannot be negative")
    if args.rigid_if_over is not None and args.rigid_if_over < 0:
        parser.error("--rigid-if-over cannot be negative")
    return args


def discover_inputs(path: Path, recursive: bool) -> tuple[list[Path], Path]:
    if path.is_file():
        if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
            raise ValueError(f"Unsupported ligand format: {path.suffix}")
        return [path], path.parent
    if not path.is_dir():
        raise FileNotFoundError(f"Input does not exist: {path}")
    iterator = path.rglob("*") if recursive else path.glob("*")
    files = sorted(
        item.resolve()
        for item in iterator
        if item.is_file() and item.suffix.lower() in SUPPORTED_EXTENSIONS
    )
    if not files:
        raise RuntimeError(f"No supported ligand files found under {path}")
    return files, path


def output_bases(files: list[Path], input_root: Path) -> dict[Path, tuple[Path, str]]:
    records: list[tuple[Path, Path, str]] = []
    counts: dict[tuple[str, str], int] = {}
    for source in files:
        relative_parent = source.parent.relative_to(input_root)
        stem = safe_stem(source.stem)
        key = (str(relative_parent).casefold(), stem.casefold())
        counts[key] = counts.get(key, 0) + 1
        records.append((source, relative_parent, stem))
    output: dict[Path, tuple[Path, str]] = {}
    for source, relative_parent, stem in records:
        key = (str(relative_parent).casefold(), stem.casefold())
        if counts[key] > 1:
            stem += "_" + safe_stem(source.suffix.lstrip("."), "file")
        output[source] = (relative_parent, stem)
    return output


def split_and_protonate(
    source: Path, staging: Path, obabel: str, ph: float
) -> list[Path]:
    split_template = staging / "molecule.sdf"
    command = [
        obabel,
        str(source),
        "-O",
        str(split_template),
        "-m",
        "-p",
        str(ph),
    ]
    completed = run_command(command, timeout=1800)
    generated = sorted(staging.glob("molecule*.sdf"))
    if completed.returncode != 0 or not generated:
        raise RuntimeError(command_error(completed) or "Open Babel produced no molecules")
    return generated


def generate_3d(
    input_sdf: Path,
    output_sdf: Path,
    *,
    seed: int,
    optimize_steps: int,
    keep_salts: bool,
    keep_existing_3d: bool,
) -> tuple[int, str]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False)
    molecule = next((item for item in supplier if item is not None), None)
    if molecule is None:
        raise RuntimeError("RDKit could not parse the protonated SDF")
    if not keep_salts:
        fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
        if fragments:
            molecule = max(fragments, key=lambda item: item.GetNumHeavyAtoms())
    molecule = Chem.AddHs(molecule, addCoords=True)
    has_3d = molecule.GetNumConformers() > 0 and molecule.GetConformer().Is3D()
    method = "existing"
    if not (keep_existing_3d and has_3d):
        molecule.RemoveAllConformers()
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = seed
        parameters.useRandomCoords = True
        status = AllChem.EmbedMolecule(molecule, parameters)
        if status != 0:
            raise RuntimeError("RDKit ETKDG could not generate a 3D conformer")
        method = "ETKDGv3"
    if optimize_steps:
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule, maxIters=optimize_steps)
            method += "+MMFF94"
        elif AllChem.UFFHasAllMoleculeParams(molecule):
            AllChem.UFFOptimizeMolecule(molecule, maxIters=optimize_steps)
            method += "+UFF"
        else:
            method += "+unoptimized"
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(molecule)
    writer.close()
    if not output_sdf.is_file() or output_sdf.stat().st_size == 0:
        raise RuntimeError("RDKit did not write the 3D SDF")
    return molecule.GetNumAtoms(), method


def rigidify_pdbqt(path: Path) -> None:
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    remarks = [line for line in lines if line.startswith("REMARK")]
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms:
        raise RuntimeError("Cannot rigidify a PDBQT that has no atoms")
    output = remarks + ["REMARK RIGIDIFIED_BY DOCKING_PREPROCESSING", "ROOT"]
    output.extend(atoms)
    output.extend(["ENDROOT", "TORSDOF 0", ""])
    temporary = path.with_name(path.stem + ".rigid.tmp.pdbqt")
    temporary.write_text("\n".join(output), encoding="ascii")
    temporary.replace(path)


def main() -> int:
    args = parse_args()
    require_molecular_docking_environment()
    obabel = find_command("obabel")
    meeko = find_command("mk_prepare_ligand")
    source_path = Path(args.input).resolve()
    output_root = Path(args.output_root).resolve()
    files, input_root = discover_inputs(source_path, not args.no_recursive)
    if args.limit is not None:
        files = files[: args.limit]
    bases = output_bases(files, input_root)
    output_root.mkdir(parents=True, exist_ok=True)
    parameters = {
        "ph": args.ph,
        "seed": args.seed,
        "optimize_steps": args.optimize_steps,
        "keep_salts": args.keep_salts,
        "keep_existing_3d": args.keep_existing_3d,
        "rigid_macrocycles": args.rigid_macrocycles,
        "rigid_if_over": args.rigid_if_over,
        "charge_model": args.charge_model,
    }
    versions = tool_versions()

    def prepare_source(source: Path) -> list[dict]:
        started = time.perf_counter()
        relative_parent, base = bases[source]
        target_dir = output_root / relative_parent
        target_dir.mkdir(parents=True, exist_ok=True)
        rows: list[dict] = []
        try:
            with tempfile.TemporaryDirectory(prefix="docking_ligand_prep_") as temp_name:
                staging = Path(temp_name)
                molecules = split_and_protonate(source, staging, obabel, args.ph)
                for index, protonated in enumerate(molecules, 1):
                    name = base if len(molecules) == 1 else f"{base}_{index:03d}"
                    target = target_dir / f"{name}.pdbqt"
                    metadata_path = target_dir / f"{name}_metadata.json"
                    row = {
                        "source": str(source),
                        "molecule_index": index,
                        "name": name,
                        "status": "failed",
                        "output_pdbqt": "",
                        "torsdof": "",
                        "atom_count": "",
                        "rigidified": False,
                        "elapsed_seconds": "",
                        "message": "",
                    }
                    try:
                        if target.exists() and not args.overwrite:
                            raise FileExistsError(
                                f"Output exists; pass --overwrite: {target}"
                            )
                        prepared_sdf = staging / f"prepared_{index:03d}.sdf"
                        atom_count, conformer_method = generate_3d(
                            protonated,
                            prepared_sdf,
                            seed=args.seed + index - 1,
                            optimize_steps=args.optimize_steps,
                            keep_salts=args.keep_salts,
                            keep_existing_3d=args.keep_existing_3d,
                        )
                        temporary_pdbqt = target_dir / f"{name}.tmp.pdbqt"
                        temporary_pdbqt.unlink(missing_ok=True)
                        command = [
                            meeko,
                            "-i",
                            str(prepared_sdf),
                            "-o",
                            str(temporary_pdbqt),
                            "--charge_model",
                            args.charge_model,
                        ]
                        if args.rigid_macrocycles:
                            command.append("--rigid_macrocycles")
                        completed = run_command(command, timeout=1800)
                        if completed.returncode != 0 or not temporary_pdbqt.is_file():
                            raise RuntimeError(
                                command_error(completed) or "Meeko produced no PDBQT"
                            )
                        if pdbqt_atom_count(temporary_pdbqt) == 0:
                            raise RuntimeError("Generated PDBQT contains no atoms")
                        torsdof = pdbqt_torsdof(temporary_pdbqt)
                        rigidified = (
                            args.rigid_if_over is not None
                            and torsdof is not None
                            and torsdof > args.rigid_if_over
                        )
                        if rigidified:
                            rigidify_pdbqt(temporary_pdbqt)
                            torsdof = 0
                        temporary_pdbqt.replace(target)
                        if args.keep_intermediate:
                            shutil.copy2(prepared_sdf, target_dir / f"{name}_prepared.sdf")
                        metadata = {
                            "source": str(source),
                            "source_sha256": sha256(source),
                            "molecule_index": index,
                            "parameters": parameters,
                            "conformer_method": conformer_method,
                            "input_atom_count_with_hydrogens": atom_count,
                            "output_pdbqt": str(target),
                            "output_sha256": sha256(target),
                            "output_atom_count": pdbqt_atom_count(target),
                            "torsdof": torsdof,
                            "rigidified": rigidified,
                            "tools": versions,
                        }
                        write_json(metadata_path, metadata)
                        row.update(
                            status="completed_rigidified" if rigidified else "completed",
                            output_pdbqt=str(target),
                            torsdof=torsdof,
                            atom_count=pdbqt_atom_count(target),
                            rigidified=rigidified,
                        )
                    except Exception as exc:
                        row["message"] = str(exc).replace("\t", " ").replace("\n", " ")
                    row["elapsed_seconds"] = f"{time.perf_counter() - started:.2f}"
                    rows.append(row)
        except (Exception, subprocess.TimeoutExpired) as exc:
            rows.append(
                {
                    "source": str(source),
                    "molecule_index": "",
                    "name": base,
                    "status": "failed",
                    "output_pdbqt": "",
                    "torsdof": "",
                    "atom_count": "",
                    "rigidified": False,
                    "elapsed_seconds": f"{time.perf_counter() - started:.2f}",
                    "message": str(exc).replace("\t", " ").replace("\n", " "),
                }
            )
        return rows

    rows: list[dict] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        future_map = {pool.submit(prepare_source, source): source for source in files}
        for completed_count, future in enumerate(
            concurrent.futures.as_completed(future_map), 1
        ):
            source = future_map[future]
            try:
                source_rows = future.result()
            except Exception as exc:
                source_rows = [{field: "" for field in FIELDS}]
                source_rows[0].update(
                    source=str(source), status="failed", message=repr(exc)
                )
            rows.extend(source_rows)
            rows.sort(key=lambda item: (item["source"], str(item["molecule_index"])))
            write_tsv(output_root / "preprocessing_results.tsv", rows, FIELDS)
            statuses = ", ".join(row["status"] for row in source_rows)
            print(f"[{completed_count}/{len(files)}] {source.name}: {statuses}", flush=True)

    failed = [row for row in rows if not str(row["status"]).startswith("completed")]
    summary = {
        "input_files": len(files),
        "output_molecules": len(rows),
        "completed": len(rows) - len(failed),
        "failed": len(failed),
        "parameters": parameters,
        "tools": versions,
        "results_tsv": str(output_root / "preprocessing_results.tsv"),
    }
    write_json(output_root / "preprocessing_summary.json", summary)
    print(f"Results: {output_root / 'preprocessing_results.tsv'}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
