from __future__ import annotations

import argparse
import concurrent.futures
import json
import subprocess
import time
from pathlib import Path

from workflow_common import (
    build_complex,
    find_command,
    find_vina,
    load_config,
    numeric_key,
    package_path,
    parse_case_list,
    pdbqt_atom_count,
    pdbqt_torsdof,
    run_command,
    sha256,
    vina_scores,
    write_tsv,
)


FIELDS = [
    "group",
    "number",
    "status",
    "input_pdbqt",
    "input_sha256",
    "best_affinity_kcal_mol",
    "torsdof",
    "pose_count",
    "elapsed_seconds",
    "all_poses_pdbqt",
    "output_sha256",
    "best_pose_pdb",
    "best_complex_pdb",
    "log",
    "message",
]


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def result_fingerprint(
    ligand: Path,
    receptor: Path,
    display_receptor: Path,
    vina: Path,
    docking: dict,
) -> dict:
    return {
        "input_sha256": sha256(ligand),
        "receptor_sha256": sha256(receptor),
        "display_receptor_sha256": sha256(display_receptor),
        "vina_sha256": sha256(vina),
        "docking": docking,
    }


def load_metadata(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def save_metadata(path: Path, metadata: dict) -> None:
    atomic_write_text(path, json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def subprocess_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description="Run deterministic batch docking with Vina.")
    parser.add_argument("--config")
    parser.add_argument("--input-root", help="Prepared ligand PDBQT root.")
    parser.add_argument("--output-root")
    parser.add_argument("--groups", help="Comma-separated groups")
    parser.add_argument(
    "--cases",
    help=(
        "Comma-separated cases such as "
        "DEFAULT/ligand_001,DEFAULT/ligand_002"
    ),
)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--allow-high-torsion",
        action="store_true",
        help="Allow ligands above max_flexible_torsions (may be very slow).",
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be at least 1")

    config = load_config(args.config)
    if args.input_root:
        input_root = Path(args.input_root).resolve()
    else:
        input_root = package_path(config["prepared_ligand_root"])
    output_root = (
        Path(args.output_root).resolve()
        if args.output_root
        else package_path(config["default_output_root"])
    )
    configured_groups = list(config["groups"])
    selected_groups = (
        {item.strip() for item in args.groups.split(",") if item.strip()}
        if args.groups
        else set(configured_groups)
    )
    unknown_groups = selected_groups - set(configured_groups)
    if unknown_groups:
        parser.error("Unknown group(s): " + ", ".join(sorted(unknown_groups)))
    selected_cases = parse_case_list(args.cases)
    invalid_cases = set(config.get("invalid_cases", []))
    unsupported_cases = set(config.get("unsupported_cases", []))
    vina = find_vina(config)
    obabel = find_command("obabel")
    receptor = package_path(config["prepared_receptor"])
    display_receptor = package_path(config["display_receptor"])
    docking = config["docking"]
    for label, path in (
        ("ligand root", input_root),
        ("prepared receptor", receptor),
        ("display receptor", display_receptor),
    ):
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")
    if not input_root.is_dir():
        raise NotADirectoryError(f"Ligand root is not a directory: {input_root}")

    jobs = []
    discovered_cases = set()
    for group in configured_groups:
        if group not in selected_groups:
            continue
        for ligand in sorted((input_root / group).glob("*.pdbqt"), key=numeric_key):
            case = f"{group}/{ligand.stem}"
            discovered_cases.add(case)
            if selected_cases is not None and case not in selected_cases:
                continue
            jobs.append((group, ligand))
    if selected_cases is not None:
        missing_cases = selected_cases - discovered_cases
        if missing_cases:
            parser.error("Unknown or unavailable case(s): " + ", ".join(sorted(missing_cases)))
    if args.limit is not None:
        jobs = jobs[: args.limit]
    if not jobs:
        raise RuntimeError(f"No PDBQT ligands selected under {input_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    results_path = output_root / "results.tsv"
    summary_path = output_root / "SUMMARY.md"
    write_tsv(results_path, [], FIELDS)
    atomic_write_text(
        summary_path,
        "# Docking run summary\n\n- Status: in progress\n"
        f"- Selected inputs: {len(jobs)}\n",
    )

    print(f"Vina: {vina}")
    print(f"Receptor: {receptor}")
    print(f"Ligand root: {input_root}")
    print(f"Output root: {output_root}")
    print(f"Jobs: {len(jobs)}, workers: {args.workers}, cpu/job: {docking['cpu_per_job']}")

    def dock_one(group: str, ligand: Path) -> dict:
        number = ligand.stem
        case = f"{group}/{number}"
        result_dir = output_root / group
        result_dir.mkdir(parents=True, exist_ok=True)
        all_poses = result_dir / f"{number}_all_poses.pdbqt"
        best_pose = result_dir / f"{number}_best_pose.pdb"
        best_complex = result_dir / f"{number}_best_complex.pdb"
        log = result_dir / f"{number}_vina.log"
        metadata_path = result_dir / f"{number}_metadata.json"
        temporary_poses = result_dir / f"{number}_all_poses.tmp.pdbqt"
        temporary_pose = result_dir / f"{number}_best_pose.tmp.pdb"
        temporary_complex = result_dir / f"{number}_best_complex.tmp.pdb"
        started = time.perf_counter()
        status = "completed"
        message = ""
        log_available = False
        input_hash = sha256(ligand)
        torsdof = pdbqt_torsdof(ligand)

        def make_row(scores: list[float], *, include_outputs: bool) -> dict:
            best = scores[0] if scores else None
            elapsed = time.perf_counter() - started
            return {
                "group": group,
                "number": number,
                "status": status,
                "input_pdbqt": str(ligand),
                "input_sha256": input_hash,
                "best_affinity_kcal_mol": "" if best is None else f"{best:.3f}",
                "torsdof": "" if torsdof is None else torsdof,
                "pose_count": len(scores),
                "elapsed_seconds": f"{elapsed:.2f}",
                "all_poses_pdbqt": str(all_poses) if (include_outputs or status == "conversion_failed") and all_poses.is_file() else "",
                "output_sha256": sha256(all_poses) if (include_outputs or status == "conversion_failed") and all_poses.is_file() else "",
                "best_pose_pdb": str(best_pose) if include_outputs and best_pose.is_file() else "",
                "best_complex_pdb": str(best_complex) if include_outputs and best_complex.is_file() else "",
                "log": str(log) if log_available and log.is_file() else "",
                "message": message.replace("\t", " ").replace("\n", " "),
            }

        if case in invalid_cases:
            status = "skipped_invalid"
            message = "Case is listed in config invalid_cases"
            return make_row([], include_outputs=False)
        if case in unsupported_cases:
            status = "skipped_unsupported"
            message = "Case is listed in config unsupported_cases"
            return make_row([], include_outputs=False)
        if ligand.stat().st_size == 0 or pdbqt_atom_count(ligand) == 0:
            status = "invalid_input"
            message = "Ligand PDBQT contains no ATOM/HETATM records"
            return make_row([], include_outputs=False)
        max_torsions = config.get("max_flexible_torsions")
        if (
            not args.allow_high_torsion
            and max_torsions is not None
            and torsdof is not None
            and torsdof > max_torsions
        ):
            status = "requires_rigid_input"
            message = (
                f"TORSDOF {torsdof} exceeds configured maximum {max_torsions}; "
                "prepare a rigid PDBQT or pass --allow-high-torsion"
            )
            return make_row([], include_outputs=False)

        fingerprint = result_fingerprint(
            ligand, receptor, display_receptor, vina, docking
        )
        metadata = load_metadata(metadata_path)
        resume_valid = (
            args.resume
            and metadata.get("fingerprint") == fingerprint
            and all_poses.is_file()
            and metadata.get("all_poses_sha256") == sha256(all_poses)
            and bool(vina_scores(all_poses))
        )
        if resume_valid:
            status = "resumed"
            log_available = log.is_file()
            derived_valid = (
                best_pose.is_file()
                and best_complex.is_file()
                and metadata.get("best_pose_sha256") == sha256(best_pose)
                and metadata.get("best_complex_sha256") == sha256(best_complex)
            )
            if derived_valid:
                return make_row(vina_scores(all_poses), include_outputs=True)
        else:
            center = docking["center"]
            size = docking["size"]
            temporary_poses.unlink(missing_ok=True)
            command = [
                str(vina),
                "--receptor", str(receptor),
                "--ligand", str(ligand),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--exhaustiveness", str(docking["exhaustiveness"]),
                "--num_modes", str(docking["num_modes"]),
                "--energy_range", str(docking["energy_range"]),
                "--seed", str(docking["seed"]),
                "--cpu", str(docking["cpu_per_job"]),
                "--out", str(temporary_poses),
            ]
            try:
                completed = run_command(command, timeout=7200)
                return_code = completed.returncode
                stdout = completed.stdout
                stderr = completed.stderr
            except subprocess.TimeoutExpired as exc:
                return_code = "TIMEOUT"
                stdout = subprocess_text(exc.stdout)
                stderr = subprocess_text(exc.stderr)
            atomic_write_text(
                log,
                "COMMAND\n" + "\n".join(command) + "\n\nSTDOUT\n" + stdout
                + "\nSTDERR\n" + stderr
                + f"\nRETURN_CODE\n{return_code}\n",
            )
            log_available = True
            if return_code != 0 or not temporary_poses.is_file():
                status = "docking_failed"
                message = stderr.strip() or f"Vina exit code {return_code}"
                temporary_poses.unlink(missing_ok=True)
                return make_row([], include_outputs=False)
            temporary_scores = vina_scores(temporary_poses)
            if not temporary_scores:
                status = "no_scores"
                message = "No valid REMARK VINA RESULT records"
                temporary_poses.unlink(missing_ok=True)
                return make_row([], include_outputs=False)
            temporary_poses.replace(all_poses)
            metadata = {
                "fingerprint": fingerprint,
                "all_poses_sha256": sha256(all_poses),
            }
            save_metadata(metadata_path, metadata)
        scores = vina_scores(all_poses) if all_poses.is_file() else []
        if status in {"completed", "resumed"} and not scores:
            status = "no_scores"
            message = "No REMARK VINA RESULT records"
            return make_row([], include_outputs=False)

        if status in {"completed", "resumed"}:
            temporary_pose.unlink(missing_ok=True)
            temporary_complex.unlink(missing_ok=True)
            conversion_command = [
                obabel,
                str(all_poses),
                "-O",
                str(temporary_pose),
                "-f",
                "1",
                "-l",
                "1",
            ]
            try:
                converted = run_command(conversion_command, timeout=600)
                conversion_return_code = converted.returncode
                conversion_error = converted.stderr
            except subprocess.TimeoutExpired as exc:
                conversion_return_code = "TIMEOUT"
                conversion_error = subprocess_text(exc.stderr)
            if conversion_return_code == 0 and temporary_pose.is_file():
                build_complex(display_receptor, temporary_pose, temporary_complex)
                temporary_pose.replace(best_pose)
                temporary_complex.replace(best_complex)
                metadata.update(
                    best_pose_sha256=sha256(best_pose),
                    best_complex_sha256=sha256(best_complex),
                )
                save_metadata(metadata_path, metadata)
            else:
                status = "conversion_failed"
                detail = conversion_error.strip() or f"Open Babel exit code {conversion_return_code}"
                message = "Docking completed; best-pose PDB conversion failed: " + detail
        return make_row(scores, include_outputs=status in {"completed", "resumed"})

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        future_map = {
            pool.submit(dock_one, group, ligand): (group, ligand)
            for group, ligand in jobs
        }
        for index, future in enumerate(concurrent.futures.as_completed(future_map), 1):
            group, ligand = future_map[future]
            try:
                row = future.result()
            except Exception as exc:
                row = {field: "" for field in FIELDS}
                row.update(
                    group=group,
                    number=ligand.stem,
                    status="exception",
                    input_pdbqt=str(ligand),
                    message=repr(exc),
                )
            rows.append(row)
            partial_rows = sorted(
                rows,
                key=lambda item: (
                    configured_groups.index(item["group"]),
                    numeric_key(Path(item["number"])),
                ),
            )
            write_tsv(results_path, partial_rows, FIELDS)
            print(
                f"[{index}/{len(jobs)}] {group}/{ligand.stem} "
                f"{row['status']} score={row['best_affinity_kcal_mol']} "
                f"torsdof={row['torsdof']}",
                flush=True,
            )

    group_order = {group: index for index, group in enumerate(configured_groups)}
    rows.sort(key=lambda row: (group_order[row["group"]], numeric_key(Path(row["number"]))))
    write_tsv(results_path, rows, FIELDS)
    success_statuses = {"completed", "resumed"}
    skipped_statuses = {"skipped_invalid", "skipped_unsupported"}
    completed_rows = [row for row in rows if row["status"] in success_statuses]
    skipped_rows = [row for row in rows if row["status"] in skipped_statuses]
    failed_rows = [
        row
        for row in rows
        if row["status"] not in success_statuses | skipped_statuses
    ]
    ranked = sorted(
        [row for row in completed_rows if row["best_affinity_kcal_mol"]],
        key=lambda row: float(row["best_affinity_kcal_mol"]),
    )
    summary = [
        "# Docking run summary",
        "",
        f"- Inputs: {len(rows)}",
        f"- Completed/resumed: {len(completed_rows)}",
        f"- Configured skips: {len(skipped_rows)}",
        f"- Failed/action required: {len(failed_rows)}",
        f"- Ligand root: `{input_root}`",
        f"- Receptor: `{receptor}`",
        f"- Vina: `{vina}`",
        "",
        "| Rank | Group | Number | Score | TORSDOF | Poses | Status |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for rank, row in enumerate(ranked, 1):
        summary.append(
            f"| {rank} | {row['group']} | {row['number']} | "
            f"{row['best_affinity_kcal_mol']} | {row['torsdof']} | "
            f"{row['pose_count']} | {row['status']} |"
        )
    non_completed = [row for row in rows if row["status"] not in success_statuses]
    if non_completed:
        summary.extend(
            [
                "",
                "## Non-completed inputs",
                "",
                "| Group | Number | TORSDOF | Status | Message |",
                "|---|---:|---:|---|---|",
            ]
        )
        for row in non_completed:
            message = row["message"].replace("|", "\\|")
            summary.append(
                f"| {row['group']} | {row['number']} | {row['torsdof']} | "
                f"{row['status']} | {message} |"
            )
    atomic_write_text(summary_path, "\n".join(summary) + "\n")
    print(f"Results: {results_path}")
    return 0 if not failed_rows else 1


if __name__ == "__main__":
    raise SystemExit(main())
