from __future__ import annotations

import argparse
import re
import shutil
import tempfile
import time
from pathlib import Path

from preprocessing_common import (
    PREPROCESSING_ROOT,
    command_error,
    find_command,
    pdbqt_atom_count,
    require_molecular_docking_environment,
    run_command,
    safe_stem,
    sha256,
    tool_versions,
    write_json,
)


SUPPORTED_EXTENSIONS = {".pdb", ".ent", ".cif", ".mmcif"}
WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD", "SOL"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare a raw receptor structure as an AutoDock Vina PDBQT."
    )
    parser.add_argument("--input", required=True, help="PDB/mmCIF receptor file")
    parser.add_argument(
        "--output-root",
        default=str(PREPROCESSING_ROOT / "output" / "receptor"),
    )
    parser.add_argument("--name", help="Output basename; defaults to the input stem")
    parser.add_argument("--chains", help="Comma-separated chain IDs to retain")
    parser.add_argument("--model", type=int, default=1)
    parser.add_argument("--altloc", default="A", help="Preferred alternate location")
    parser.add_argument("--keep-hetero", action="store_true")
    parser.add_argument(
        "--keep-residue",
        action="append",
        default=[],
        metavar="RESNAME",
        help="HETATM residue name to retain; may be repeated or comma-separated",
    )
    parser.add_argument("--keep-water", action="store_true")
    parser.add_argument("--keep-input-hydrogens", action="store_true")
    parser.add_argument("--allow-bad-res", action="store_true")
    parser.add_argument(
        "--set-template",
        action="append",
        default=[],
        metavar="SPEC",
        help="Meeko residue template assignment, e.g. A:42=HID",
    )
    parser.add_argument("--delete-residues", metavar="SPEC")
    parser.add_argument("--add-templates", metavar="JSON_OR_RESNAME:SDF")
    parser.add_argument(
        "--charge-model",
        choices=("gasteiger", "zero"),
        default="gasteiger",
    )
    parser.add_argument("--box-center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--box-size", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--box-enveloping", help="Molecule whose atoms define the box")
    parser.add_argument("--padding", type=float, default=5.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-intermediate", action="store_true")
    args = parser.parse_args()
    if args.model < 1:
        parser.error("--model must be at least 1")
    if (args.box_center is None) != (args.box_size is None):
        parser.error("--box-center and --box-size must be supplied together")
    if args.box_enveloping and args.box_center:
        parser.error("Use either an explicit box or --box-enveloping, not both")
    return args


def atom_is_hydrogen(line: str) -> bool:
    padded = line.ljust(80)
    element = padded[76:78].strip().upper()
    if element:
        return element in {"H", "D"}
    atom_name = re.sub(r"^[0-9]+", "", padded[12:16].strip()).upper()
    return atom_name.startswith(("H", "D"))


def occupancy(line: str) -> float:
    try:
        return float(line.ljust(60)[54:60])
    except ValueError:
        return 0.0


def clean_receptor_pdb(
    source: Path,
    output: Path,
    *,
    chains: set[str] | None,
    model_number: int,
    preferred_altloc: str,
    keep_hetero: bool,
    keep_residues: set[str],
    keep_water: bool,
    keep_hydrogens: bool,
) -> dict[str, int]:
    candidates: dict[tuple[str, str, str, str, str], list[str]] = {}
    residue_order: dict[tuple[str, str, str, str], int] = {}
    atom_order: dict[tuple[str, str, str, str, str], int] = {}
    stats = {
        "atom_records_seen": 0,
        "atoms_retained": 0,
        "dropped_chain": 0,
        "dropped_water": 0,
        "dropped_hetero": 0,
        "dropped_hydrogen": 0,
        "alternate_locations_removed": 0,
    }
    has_models = False
    active_model = True
    for raw_line in source.read_text(encoding="ascii", errors="replace").splitlines():
        if raw_line.startswith("MODEL"):
            has_models = True
            try:
                active_model = int(raw_line.split()[1]) == model_number
            except (IndexError, ValueError):
                active_model = model_number == 1
            continue
        if raw_line.startswith("ENDMDL"):
            if has_models and active_model:
                break
            active_model = False
            continue
        if not active_model or not raw_line.startswith(("ATOM  ", "HETATM")):
            continue
        stats["atom_records_seen"] += 1
        line = raw_line.ljust(80)
        record = line[:6]
        atom_name = line[12:16]
        residue = line[17:20].strip().upper()
        chain = line[21].strip()
        residue_number = line[22:26]
        insertion_code = line[26]
        if chains is not None and chain not in chains:
            stats["dropped_chain"] += 1
            continue
        if residue in WATER_RESIDUES and not keep_water:
            stats["dropped_water"] += 1
            continue
        if record == "HETATM" and not (
            keep_hetero or residue in keep_residues or (keep_water and residue in WATER_RESIDUES)
        ):
            stats["dropped_hetero"] += 1
            continue
        if not keep_hydrogens and atom_is_hydrogen(line):
            stats["dropped_hydrogen"] += 1
            continue
        key = (chain, residue_number, insertion_code, residue, atom_name)
        residue_key = key[:4]
        if residue_key not in residue_order:
            residue_order[residue_key] = len(residue_order)
        if key not in atom_order:
            atom_order[key] = len(atom_order)
        candidates.setdefault(key, []).append(line)

    chosen_items: list[tuple[int, int, str]] = []
    preferred = preferred_altloc[:1].upper() if preferred_altloc else "A"
    for key, lines in candidates.items():
        selected = max(
            lines,
            key=lambda line: (
                line[16].upper() == preferred,
                line[16] == " ",
                occupancy(line),
            ),
        )
        stats["alternate_locations_removed"] += len(lines) - 1
        selected = selected[:16] + " " + selected[17:]
        residue_key = key[:4]
        chosen_items.append(
            (residue_order[residue_key], atom_order[key], selected.rstrip())
        )
    chosen_items.sort(key=lambda item: (item[0], item[1]))
    chosen = [item[2] for item in chosen_items]
    if not chosen:
        raise RuntimeError("No receptor atoms remain after cleaning")

    output_lines: list[str] = []
    last_chain: str | None = None
    for line in chosen:
        current_chain = line.ljust(22)[21]
        if last_chain is not None and current_chain != last_chain:
            output_lines.append("TER")
        output_lines.append(line)
        last_chain = current_chain
    output_lines.extend(["TER", "END", ""])
    output.write_text("\n".join(output_lines), encoding="ascii")
    stats["atoms_retained"] = len(chosen)
    return stats


def mmcif_to_pdb(source: Path, output: Path) -> None:
    from Bio.PDB import MMCIFParser, PDBIO

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("receptor", str(source))
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(output))
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("Biopython did not write a PDB from the mmCIF input")


def add_template_polar_hydrogens(source: Path, output: Path) -> int:
    """Add canonical Ser/Thr/Tyr/Cys side-chain hydrogens before Meeko bond perception."""
    from rdkit import Chem

    protein_only = output.with_name(output.stem + ".protein_only.pdb")
    protein_lines = [
        line
        for line in source.read_text(encoding="ascii", errors="replace").splitlines()
        if line.startswith("ATOM  ")
    ]
    protein_only.write_text(
        "\n".join(protein_lines + ["TER", "END", ""]), encoding="ascii"
    )
    molecule = Chem.MolFromPDBFile(
        str(protein_only), removeHs=False, sanitize=False, proximityBonding=True
    )
    protein_only.unlink(missing_ok=True)
    if molecule is None:
        raise RuntimeError("RDKit could not parse the cleaned receptor PDB")
    try:
        Chem.SanitizeMol(molecule)
        heavy_atom_count = molecule.GetNumAtoms()
        hydrogenated = Chem.AddHs(
            molecule, addCoords=True, addResidueInfo=True
        )
    except Exception as exc:
        raise RuntimeError(
            "RDKit could not place receptor polar hydrogens; provide a "
            "hydrogenated PDB and pass --keep-input-hydrogens"
        ) from exc

    canonical = {
        ("SER", "OG"): "HG",
        ("THR", "OG1"): "HG1",
        ("TYR", "OH"): "HH",
        ("CYS", "SG"): "HG",
    }
    conformer = hydrogenated.GetConformer()
    additions: list[str] = []
    serial = 90000
    for atom in hydrogenated.GetAtoms():
        if atom.GetIdx() < heavy_atom_count or atom.GetAtomicNum() != 1:
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 1:
            continue
        parent_info = neighbors[0].GetPDBResidueInfo()
        if parent_info is None:
            continue
        residue_name = parent_info.GetResidueName().strip().upper()
        parent_name = parent_info.GetName().strip().upper()
        hydrogen_name = canonical.get((residue_name, parent_name))
        if hydrogen_name is None:
            continue
        position = conformer.GetAtomPosition(atom.GetIdx())
        chain = (parent_info.GetChainId() or " ")[:1]
        insertion_code = (parent_info.GetInsertionCode() or " ")[:1]
        residue_number = parent_info.GetResidueNumber()
        additions.append(
            f"ATOM  {serial:5d} {hydrogen_name:>4s} {residue_name:>3s} "
            f"{chain}{residue_number:4d}{insertion_code}   "
            f"{position.x:8.3f}{position.y:8.3f}{position.z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {'H':>2s}  "
        )
        serial += 1

    unordered = output.with_name(output.stem + ".unordered.pdb")
    atom_lines = [
        line
        for line in source.read_text(encoding="ascii", errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    unordered.write_text(
        "\n".join(atom_lines + additions + ["TER", "END", ""]), encoding="ascii"
    )
    clean_receptor_pdb(
        unordered,
        output,
        chains=None,
        model_number=1,
        preferred_altloc="A",
        keep_hetero=True,
        keep_residues=set(),
        keep_water=True,
        keep_hydrogens=True,
    )
    unordered.unlink(missing_ok=True)
    return len(additions)


def main() -> int:
    args = parse_args()
    require_molecular_docking_environment()
    source = Path(args.input).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Receptor input does not exist: {source}")
    if source.suffix.lower() not in SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported receptor format: {source.suffix}")
    output_root = Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    name = safe_stem(args.name or source.stem, "receptor")
    output_pdbqt = output_root / f"{name}_prep.pdbqt"
    output_pdb = output_root / f"{name}_prep.pdb"
    metadata_path = output_root / f"{name}_metadata.json"
    box_path = output_root / f"{name}_vina_box.txt"
    expected = [output_pdbqt, output_pdb, metadata_path]
    if args.box_center or args.box_enveloping:
        expected.append(box_path)
    existing = [path for path in expected if path.exists()]
    if existing and not args.overwrite:
        joined = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Output exists; pass --overwrite: {joined}")

    chains = None
    if args.chains:
        chains = {item.strip() for item in args.chains.split(",")}
    keep_residues = {
        item.strip().upper()
        for value in args.keep_residue
        for item in value.split(",")
        if item.strip()
    }
    meeko = find_command("mk_prepare_receptor")
    started = time.perf_counter()

    with tempfile.TemporaryDirectory(prefix="docking_receptor_prep_") as temp_name:
        staging = Path(temp_name)
        source_pdb = staging / "source.pdb"
        if source.suffix.lower() in {".pdb", ".ent"}:
            shutil.copy2(source, source_pdb)
        else:
            mmcif_to_pdb(source, source_pdb)

        cleaned = staging / "cleaned.pdb"
        cleaning_stats = clean_receptor_pdb(
            source_pdb,
            cleaned,
            chains=chains,
            model_number=args.model,
            preferred_altloc=args.altloc,
            keep_hetero=args.keep_hetero,
            keep_residues=keep_residues,
            keep_water=args.keep_water,
            keep_hydrogens=args.keep_input_hydrogens,
        )
        receptor_for_meeko = cleaned
        polar_hydrogens_added = 0
        if not args.keep_input_hydrogens:
            receptor_for_meeko = staging / "polar_hydrogenated.pdb"
            polar_hydrogens_added = add_template_polar_hydrogens(
                cleaned, receptor_for_meeko
            )
        staged_pdbqt = staging / f"{name}_prep.pdbqt"
        staged_pdb = staging / f"{name}_prep.pdb"
        command = [
            meeko,
            "--read_pdb",
            str(receptor_for_meeko),
            "-o",
            str(staging / name),
            "-p",
            str(staged_pdbqt),
            "--write_pdb",
            str(staged_pdb),
            "--charge_model",
            args.charge_model,
        ]
        if args.allow_bad_res:
            command.append("-a")
        for assignment in args.set_template:
            command.extend(["-n", assignment])
        if args.delete_residues:
            command.extend(["-d", args.delete_residues])
        if args.add_templates:
            command.extend(["--add_templates", args.add_templates])
        staged_box = staging / f"{name}_vina_box.txt"
        if args.box_center:
            command.extend(["-v", str(staged_box), "--box_center"])
            command.extend(str(value) for value in args.box_center)
            command.append("--box_size")
            command.extend(str(value) for value in args.box_size)
        elif args.box_enveloping:
            enveloping = Path(args.box_enveloping).resolve()
            if not enveloping.is_file():
                raise FileNotFoundError(f"Box-enveloping molecule not found: {enveloping}")
            command.extend(
                [
                    "-v",
                    str(staged_box),
                    "--box_enveloping",
                    str(enveloping),
                    "--padding",
                    str(args.padding),
                ]
            )
        prepared = run_command(command, timeout=1800)
        if prepared.returncode != 0 or not staged_pdbqt.is_file() or not staged_pdb.is_file():
            raise RuntimeError(command_error(prepared) or "Meeko receptor preparation failed")
        atom_count = pdbqt_atom_count(staged_pdbqt)
        if atom_count == 0:
            raise RuntimeError("Generated receptor PDBQT contains no atoms")

        shutil.copy2(staged_pdbqt, output_pdbqt)
        shutil.copy2(staged_pdb, output_pdb)
        if staged_box.exists():
            shutil.copy2(staged_box, box_path)
        if args.keep_intermediate:
            shutil.copy2(cleaned, output_root / f"{name}_cleaned.pdb")
            if receptor_for_meeko != cleaned:
                shutil.copy2(
                    receptor_for_meeko,
                    output_root / f"{name}_polar_hydrogenated.pdb",
                )

    metadata = {
        "source": str(source),
        "source_sha256": sha256(source),
        "parameters": {
            "chains": sorted(chains) if chains is not None else None,
            "model": args.model,
            "preferred_altloc": args.altloc,
            "keep_hetero": args.keep_hetero,
            "keep_residues": sorted(keep_residues),
            "keep_water": args.keep_water,
            "keep_input_hydrogens": args.keep_input_hydrogens,
            "allow_bad_res": args.allow_bad_res,
            "set_template": args.set_template,
            "delete_residues": args.delete_residues,
            "charge_model": args.charge_model,
            "box_center": args.box_center,
            "box_size": args.box_size,
            "box_enveloping": args.box_enveloping,
            "padding": args.padding,
        },
        "cleaning": {
            **cleaning_stats,
            "template_polar_hydrogens_added": polar_hydrogens_added,
        },
        "outputs": {
            "prepared_receptor": str(output_pdbqt),
            "prepared_receptor_sha256": sha256(output_pdbqt),
            "display_receptor": str(output_pdb),
            "display_receptor_sha256": sha256(output_pdb),
            "vina_box": str(box_path) if box_path.exists() else None,
            "pdbqt_atom_count": atom_count,
        },
        "tools": tool_versions(),
        "elapsed_seconds": round(time.perf_counter() - started, 2),
    }
    write_json(metadata_path, metadata)
    write_json(
        output_root / f"{name}_docking_config_snippet.json",
        {
            "prepared_receptor": str(output_pdbqt),
            "display_receptor": str(output_pdb),
        },
    )
    print(f"Prepared receptor: {output_pdbqt}")
    print(f"Display receptor: {output_pdb}")
    print(f"Metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
