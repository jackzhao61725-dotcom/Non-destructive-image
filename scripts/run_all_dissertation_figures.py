"""Run approved dissertation figure and plot generation scripts."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "results" / "reproducibility_manifest.json"
NOTEBOOK_CONFIG = "configs/notebook_v1_defaults.json"
FARADAY_RESULTS_CONFIG = "configs/dissertation_results_v1.json"
DETUNING_PLOT_CONFIG = "configs/dissertation_plots_v1.json"


@dataclass(frozen=True)
class RunStep:
    name: str
    script: str
    args: tuple[str, ...]
    configs: tuple[str, ...]
    expected_outputs: tuple[str, ...]
    result_type: str

    def command(self) -> list[str]:
        return [sys.executable, self.script, *self.args]


RUN_STEPS: tuple[RunStep, ...] = (
    RunStep(
        name="validate_notebook_sections",
        script="scripts/validate_notebook_sections.py",
        args=(),
        configs=(),
        expected_outputs=(),
        result_type="validation",
    ),
    RunStep(
        name="recover_condensate_stage",
        script="scripts/recover_notebook_condensate_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/condensate_stage/condensate_density_stage.svg",
            "results/notebook_aligned_recovery/condensate_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_phase_stage",
        script="scripts/recover_notebook_phase_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/phase_stage/scalar_phase_stage.svg",
            "results/notebook_aligned_recovery/phase_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_pci_stage",
        script="scripts/recover_notebook_pci_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/pci_stage/pci_image_stage.svg",
            "results/notebook_aligned_recovery/pci_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_dgi_stage",
        script="scripts/recover_notebook_dgi_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/dgi_stage/dgi_image_stage.svg",
            "results/notebook_aligned_recovery/dgi_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_faraday_stage",
        script="scripts/recover_notebook_faraday_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/faraday_stage/faraday_dark_field_stage.svg",
            "results/notebook_aligned_recovery/faraday_stage/faraday_dual_port_signal_stage.svg",
            "results/notebook_aligned_recovery/faraday_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_deterministic_camera_stage",
        script="scripts/recover_notebook_camera_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/camera_stage/camera_deterministic_stage.svg",
            "results/notebook_aligned_recovery/camera_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_noisy_camera_stage",
        script="scripts/recover_notebook_noisy_camera_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/noisy_camera_stage/noisy_camera_stage.svg",
            "results/notebook_aligned_recovery/noisy_camera_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_multishot_stage",
        script="scripts/recover_notebook_multishot_stage.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/multishot_stage/multishot_sequence_stage.svg",
            "results/notebook_aligned_recovery/multishot_stage/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="recover_noisy_multishot_filmstrip",
        script="scripts/recover_notebook_noisy_multishot_filmstrip.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/noisy_multishot_filmstrip/noisy_multishot_pci_filmstrip.svg",
            "results/notebook_aligned_recovery/noisy_multishot_filmstrip/metadata.json",
        ),
        result_type="notebook-aligned recovery",
    ),
    RunStep(
        name="generate_condensate_three_view",
        script="scripts/generate_condensate_three_view.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG,),
        expected_outputs=(
            "results/notebook_aligned_recovery/condensate_three_view/condensate_three_view.svg",
            "results/notebook_aligned_recovery/condensate_three_view/metadata.json",
        ),
        result_type="notebook-aligned model extension",
    ),
    RunStep(
        name="generate_faraday_optimisation_results",
        script="scripts/generate_dissertation_results.py",
        args=("--config", FARADAY_RESULTS_CONFIG),
        configs=(FARADAY_RESULTS_CONFIG,),
        expected_outputs=(
            "results/faraday_optimisation_v1/detuning_tradeoff.svg",
            "results/faraday_optimisation_v1/intensity_tradeoff.svg",
            "results/faraday_optimisation_v1/exposure_time_tradeoff.svg",
            "results/faraday_optimisation_v1/metadata.json",
        ),
        result_type="representative V1 plot",
    ),
    RunStep(
        name="generate_detuning_tradeoff_plot",
        script="scripts/generate_detuning_tradeoff_plot.py",
        args=("--config", DETUNING_PLOT_CONFIG),
        configs=(NOTEBOOK_CONFIG, DETUNING_PLOT_CONFIG),
        expected_outputs=(
            "results/dissertation_plots_v1/detuning_tradeoff/detuning_tradeoff.svg",
            "results/dissertation_plots_v1/detuning_tradeoff/detuning_tradeoff_data.csv",
            "results/dissertation_plots_v1/detuning_tradeoff/metadata.json",
        ),
        result_type="dissertation physics plot",
    ),
    RunStep(
        name="audit_linear_approximation_validity",
        script="scripts/audit_linear_approximation_validity.py",
        args=("--config", NOTEBOOK_CONFIG),
        configs=(NOTEBOOK_CONFIG, FARADAY_RESULTS_CONFIG, DETUNING_PLOT_CONFIG),
        expected_outputs=(
            "docs/linear_approximation_validity_audit.md",
            "results/linear_approximation_audit/linear_approximation_summary.json",
            "results/linear_approximation_audit/metadata.json",
        ),
        result_type="numerical validity audit",
    ),
)


PENDING_ITEMS: tuple[dict[str, str], ...] = ()


def git_commit() -> str:
    commands = [
        ["git", "rev-parse", "HEAD"],
        [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"],
    ]
    for command in commands:
        try:
            return subprocess.check_output(
                command,
                cwd=REPO_ROOT,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def run_step(step: RunStep, env: dict[str, str]) -> dict[str, Any]:
    command = step.command()
    print(f"Running {step.name}: {' '.join(command)}")
    subprocess.run(command, cwd=REPO_ROOT, env=env, check=True)
    return {
        "name": step.name,
        "script": step.script,
        "command": command,
        "configs": list(step.configs),
        "expected_outputs": list(step.expected_outputs),
        "result_type": step.result_type,
        "status": "passed",
    }


def build_manifest(completed_steps: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit_hash": git_commit(),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "repository_root": str(REPO_ROOT),
        "steps": completed_steps,
        "pending_items": list(PENDING_ITEMS),
        "parameter_provenance_rule": "No figure without parameter provenance.",
    }


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help="Path for the reproducibility manifest JSON.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned steps without running scripts or writing a manifest.",
    )
    args = parser.parse_args()

    if args.dry_run:
        for step in RUN_STEPS:
            print(f"{step.name}: {' '.join(step.command())}")
        print("Pending items:")
        for item in PENDING_ITEMS:
            print(f"- {item['name']}: {item['status']}")
        return

    env = os.environ.copy()
    env["PYTHONPATH"] = "src;."
    env["PYTHONUTF8"] = "1"

    completed_steps = [run_step(step, env) for step in RUN_STEPS]
    manifest = build_manifest(completed_steps)
    write_manifest(args.manifest, manifest)
    print(f"Wrote reproducibility manifest: {args.manifest}")


if __name__ == "__main__":
    main()
