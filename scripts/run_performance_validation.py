"""Run the maintained canonical performance-validation gate.

The gate reuses the dissertation full-multishot outputs and verifies their
parameter contract, provenance and Faraday calibration boundary.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Any



REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "performance_validation_v1.json"
SOURCE_DATA = (
    REPO_ROOT
    / "results"
    / "dissertation_plots_v1"
    / "full_multishot_accumulated_snr"
    / "full_multishot_accumulated_snr_data.csv"
)
SOURCE_SUMMARY = SOURCE_DATA.with_name("full_multishot_accumulated_snr_summary.json")
SOURCE_METADATA = SOURCE_DATA.with_name("metadata.json")
SOURCE_FARADAY_LEDGER = SOURCE_DATA.with_name("faraday_canonical_reference_at_1p5GHz.csv")
SOURCE_SVG = SOURCE_DATA.with_name("full_multishot_accumulated_snr.svg")
DYNAMIC_METADATA_PATHS = (
    REPO_ROOT / "results" / "accumulated_snr_full_physics_audit" / "metadata.json",
    SOURCE_METADATA,
)
REQUIRED_LEDGER_COLUMNS = (
    "quantity",
    "value",
    "Delta/2pi",
    "P",
    "tau",
    "imaging axis",
    "normalisation",
    "N_max model",
    "QE/read",
    "ROI",
    "kappa_F",
    "calibration status",
    "config path",
    "output path",
    "git commit",
)


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(REQUIRED_LEDGER_COLUMNS)
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _git_value(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL
    ).strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _relative_path(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def _run(command: list[str], environment: dict[str, str]) -> dict[str, Any]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return {
        "command": subprocess.list2cmdline(command),
        "returncode": completed.returncode,
        "elapsed_seconds": time.perf_counter() - started,
        "output": completed.stdout.rstrip(),
        "passed": completed.returncode == 0,
    }


def run_prerequisites(source_config: str) -> list[dict[str, Any]]:
    """Run the canonical contract commands and stop at the first failure."""

    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(("src", "."))
    environment["PYTHONUTF8"] = "1"
    commands = (
        [sys.executable, "-m", "pytest", "-q"],
        [sys.executable, "scripts/validate_notebook_sections.py"],
        [
            sys.executable,
            "scripts/generate_full_multishot_accumulated_snr.py",
            "--config",
            source_config,
        ],
    )
    snapshots = {
        path: path.read_bytes() for path in DYNAMIC_METADATA_PATHS if path.exists()
    }
    results: list[dict[str, Any]] = []
    try:
        for command in commands:
            result = _run(command, environment)
            results.append(result)
            if not result["passed"]:
                break
    finally:
        for path, content in snapshots.items():
            path.write_bytes(content)
    return results


def _source_rows(detuning_ghz: float) -> list[dict[str, str]]:
    with SOURCE_DATA.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [row for row in rows if math.isclose(float(row["detuning_ghz"]), detuning_ghz)]


def _base_ledger_row(
    contract: dict[str, Any], config_path: Path, output_path: Path, commit: str
) -> dict[str, Any]:
    return {
        "Delta/2pi": f"{contract['detuning_ghz']} GHz",
        "P": f"{contract['probe_power_mw']} mW",
        "tau": f"{contract['exposure_time_us']} us",
        "imaging axis": contract["imaging_axis"],
        "normalisation": "fixed-ROI diagonal-covariance matched-template SNR",
        "N_max model": (
            "heating + initial-density reabsorption; strict integer accepted-frame "
            "count; threshold-crossing pulse excluded"
        ),
        "QE/read": (
            f"QE={contract['quantum_efficiency']}; read={contract['read_noise_electrons_per_port_pixel']} "
            "e- rms per port camera pixel"
        ),
        "ROI": f"fixed matched ROI; {contract['matched_roi_pixels']} binned camera pixels",
        "kappa_F": "1.0 for Faraday only; placeholder",
        "calibration status": "Version 1 representative and uncalibrated",
        "config path": _relative_path(config_path),
        "output path": _relative_path(output_path),
        "git commit": commit,
    }


def evaluate_canonical_gate(
    config: dict[str, Any], config_path: Path, output_path: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Evaluate exact canonical SNR and sequence invariants from maintained outputs."""

    contract = config["canonical_gate"]
    expected = contract["expected_accumulated_matched_roi_snr"]
    rtol = float(contract["relative_tolerance"])
    atol = float(contract["absolute_tolerance"])
    rows = _source_rows(float(contract["detuning_ghz"]))
    commit = _git_value("rev-parse", "HEAD")
    base = _base_ledger_row(contract, config_path, output_path, commit)
    ledger: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []

    indexed = {(row["mode"], row["noise_model"]): row for row in rows}
    for mode, noise_values in expected.items():
        for noise_model, expected_value in noise_values.items():
            source = indexed.get((mode, noise_model))
            actual = float(source["snr_total_full"]) if source is not None else math.nan
            passed = source is not None and math.isclose(actual, expected_value, rel_tol=rtol, abs_tol=atol)
            faraday = mode.startswith("Faraday")
            calibration = source["calibration_status"] if source is not None else "missing source row"
            if faraday:
                passed = passed and "uncalibrated" in calibration and source["kappa_f"] == "1.0"
            difference = actual - expected_value
            ledger.append(
                {
                    **base,
                    "quantity": f"{mode} accumulated matched-ROI SNR ({noise_model})",
                    "value": actual,
                    "mode": mode,
                    "noise model": noise_model,
                    "expected value": expected_value,
                    "absolute difference": abs(difference),
                    "relative difference": abs(difference) / abs(expected_value),
                    "status": "PASS" if passed else "FAIL",
                    "kappa_F": "1.0 placeholder" if faraday else "not applicable",
                    "calibration status": calibration,
                }
            )
            checks.append(
                {
                    "check": f"snr::{mode}::{noise_model}",
                    "actual": actual,
                    "expected": expected_value,
                    "passed": bool(passed),
                }
            )

    representative = indexed.get(("PCI", "shot_noise_only"))
    invariants = {
        "accepted_frames": (
            int(representative["n_frames_full"]) if representative is not None else None
        ),
        "post_sequence_loss": (
            float(representative["post_sequence_depletion_fraction"])
            if representative is not None
            else None
        ),
        "next_pulse_loss": (
            float(representative["next_pulse_depletion_fraction"])
            if representative is not None
            else None
        ),
        "threshold_crossing_state_included": (
            representative["threshold_crossing_state_included"] if representative is not None else None
        ),
    }
    invariant_expectations = {
        "accepted_frames": contract["accepted_frames"],
        "post_sequence_loss": contract["post_sequence_loss"],
        "next_pulse_loss": contract["next_pulse_loss"],
        "threshold_crossing_state_included": "False",
    }
    for name, expected_value in invariant_expectations.items():
        actual = invariants[name]
        passed = (
            math.isclose(float(actual), float(expected_value), rel_tol=rtol, abs_tol=atol)
            if name in {"post_sequence_loss", "next_pulse_loss"} and actual is not None
            else actual == expected_value
        )
        checks.append(
            {
                "check": f"sequence::{name}",
                "actual": actual,
                "expected": expected_value,
                "passed": bool(passed),
            }
        )

    with SOURCE_FARADAY_LEDGER.open(encoding="utf-8", newline="") as handle:
        faraday_ledger = list(csv.DictReader(handle))
    provenance_checks = {
        "source_row_count": len(rows) == 8,
        "faraday_ledger_row_count": len(faraday_ledger) == 4,
        "faraday_rows_explicitly_uncalibrated": bool(faraday_ledger)
        and all("uncalibrated" in row["calibration_status"] for row in faraday_ledger),
        "faraday_kappa_is_placeholder_one": bool(faraday_ledger)
        and all(float(row["kappa_F"]) == 1.0 for row in faraday_ledger),
        "same_fixed_roi": bool(faraday_ledger)
        and all(int(row["roi_pixel_count"]) == contract["matched_roi_pixels"] for row in faraday_ledger),
    }
    checks.extend(
        {"check": f"provenance::{name}", "actual": value, "expected": True, "passed": value}
        for name, value in provenance_checks.items()
    )
    return ledger, checks, invariants


def _write_report(
    path: Path,
    branch: str,
    commit: str,
    commands: list[dict[str, Any]],
    checks: list[dict[str, Any]],
    rows: list[dict[str, Any]],
    invariants: dict[str, Any],
    passed: bool,
    contract: dict[str, Any],
) -> None:
    command_lines = [
        f"- `{item['command']}`: {'PASS' if item['passed'] else 'FAIL'} "
        f"({item['elapsed_seconds']:.3f} s)"
        for item in commands
    ] or ["- Prerequisites intentionally skipped; numerical contract only was evaluated."]
    table_lines = [
        "| Mode | Noise | Actual | Expected | Absolute difference |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    table_lines.extend(
        f"| {row['mode']} | {row['noise model']} | {row['value']:.15g} | "
        f"{row['expected value']:.15g} | {row['absolute difference']:.3g} |"
        for row in rows
    )
    failed = [item for item in checks if not item["passed"]]
    issue_lines = [f"- `{item['check']}` failed: actual={item['actual']}, expected={item['expected']}" for item in failed]
    if not issue_lines:
        issue_lines = ["- None at the canonical gate."]
    content = "\n".join(
        [
            "# Performance Validation V1",
            "",
            "## Current stage",
            "",
            f"Canonical gate: **{'PASS' if passed else 'FAIL'}**.",
            f"Branch: `{branch}`. Source commit: `{commit}`.",
            "",
            f"This stage verifies the maintained fixed {contract['matched_roi_pixels']}-pixel matched-ROI, evolving-cloud, heating + initial-density reabsorption sequence. Peak-pixel, resolution-element, and other estimators are not mixed into this table.",
            f"The camera model is `{contract['camera_model']}`. QE={contract['quantum_efficiency']:.7g} and read noise {contract['read_noise_electrons_per_port_pixel']:.6g} e- rms per pixel and readout are provisional dissertation screening values; neither is a measurement of the installed camera.",
            "",
            "## Validation commands",
            "",
            *command_lines,
            "",
            "## Canonical four-mode results",
            "",
            *table_lines,
            "",
            f"Accepted frames: `{invariants['accepted_frames']}` (indices `0..{invariants['accepted_frames'] - 1}`). The threshold-crossing pulse is excluded.",
            f"Post-sequence loss: `{invariants['post_sequence_loss']}`.",
            f"Next-pulse loss: `{invariants['next_pulse_loss']}`.",
            "",
            "## Calibration boundary",
            "",
            "PCI and DGI rows are Version 1 notebook-aligned numerical results. Faraday rows use the phenomenological placeholder `kappa_F=1.0`; they are uncalibrated structural comparisons, not calibrated absolute or experimental predictions.",
            "",
            "## Issue register",
            "",
            *issue_lines,
            "",
            "## Evidence",
            "",
            "Machine-readable evidence is in `results/performance_validation_v1/canonical_gate.csv` and `canonical_gate.json`. The source full-multishot CSV, Faraday ledger, metadata, hashes, command output, and exact tolerances are recorded there.",
            "",
            "Later scientific, statistical, numerical, computational, and software validation stages remain separate; this report does not collapse them into a single performance score.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def run(config_path: Path, run_commands: bool = True) -> dict[str, Any]:
    config = _load_json(config_path)
    output_dir = REPO_ROOT / config["output_directory"]
    csv_path = output_dir / "canonical_gate.csv"
    json_path = output_dir / "canonical_gate.json"
    report_path = REPO_ROOT / config["report_path"]
    commands = run_prerequisites(config["source_config"]) if run_commands else []
    prerequisites_passed = all(item["passed"] for item in commands)
    if run_commands and len(commands) != 3:
        prerequisites_passed = False

    rows, checks, invariants = evaluate_canonical_gate(config, config_path, csv_path)
    numeric_passed = all(item["passed"] for item in checks)
    passed = prerequisites_passed and numeric_passed if run_commands else numeric_passed
    branch = _git_value("branch", "--show-current")
    commit = _git_value("rev-parse", "HEAD")
    _write_csv(csv_path, rows)
    evidence = {
        "stage": "canonical_gate",
        "passed": passed,
        "branch": branch,
        "git_commit": commit,
        "config_path": _relative_path(config_path),
        "source_paths": {
            "data": _relative_path(SOURCE_DATA),
            "summary": _relative_path(SOURCE_SUMMARY),
            "metadata": _relative_path(SOURCE_METADATA),
            "faraday_ledger": _relative_path(SOURCE_FARADAY_LEDGER),
        },
        "source_sha256": {
            _relative_path(path): _sha256(path)
            for path in (SOURCE_DATA, SOURCE_SUMMARY, SOURCE_FARADAY_LEDGER, SOURCE_SVG)
        },
        "commands": commands,
        "checks": checks,
        "sequence_invariants": invariants,
        "faraday_boundary": {
            "kappa_F": 1.0,
            "calibrated_absolute_prediction": False,
            "status": config["canonical_gate"]["faraday_calibration_status"],
        },
        "dynamic_metadata_policy": (
            "Known branch/commit-only metadata files are snapshotted and restored around prerequisite runs."
        ),
        "protected_paths": config["protected_paths"],
    }
    _write_json(json_path, evidence)
    _write_report(
        report_path,
        branch,
        commit,
        commands,
        checks,
        rows,
        invariants,
        passed,
        config["canonical_gate"],
    )
    if not passed:
        raise RuntimeError("Canonical performance-validation gate failed; see canonical_gate.json")
    return {"csv": csv_path, "json": json_path, "report": report_path}

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--skip-prerequisites",
        action="store_true",
        help="Evaluate maintained outputs without rerunning pytest, notebook validation, or the generator.",
    )
    args = parser.parse_args()
    config_path = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    outputs = run(config_path, run_commands=not args.skip_prerequisites)
    for label, path in outputs.items():
        print(f"- {label}: {_relative_path(path)}")


if __name__ == "__main__":
    main()
