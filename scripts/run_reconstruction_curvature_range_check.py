"""Run one declared reconstruction curvature-weight range check."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies.curvature_range import (
    run_curvature_range_check,
    write_curvature_range_check_run,
)
from non_destructive_image.reconstruction.studies.io import load_json
from non_destructive_image.reconstruction.studies.io import file_sha256
from non_destructive_image.reconstruction.studies.provenance import (
    capture_reconstruction_provenance,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
ACTIVE_CONFIGS = (
    REPO_ROOT
    / "configs"
    / "reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port.json",
    REPO_ROOT
    / "configs"
    / "reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field.json",
)


def run(config_path: Path) -> dict[str, Path]:
    """Execute and serialize one declared range check."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    generation_provenance = capture_reconstruction_provenance(
        REPO_ROOT,
        entry_points=(Path("scripts/run_reconstruction_curvature_range_check.py"),),
    )
    check_config_sha256 = file_sha256(resolved)
    check_config = load_json(resolved)
    if file_sha256(resolved) != check_config_sha256:
        raise RuntimeError("range-check config changed while it was being loaded")
    source = Path(str(check_config["source_benchmark_config"]))
    source_path = source if source.is_absolute() else REPO_ROOT / source
    benchmark_config_sha256 = file_sha256(source_path)
    benchmark_config = load_json(source_path)
    if file_sha256(source_path) != benchmark_config_sha256:
        raise RuntimeError("benchmark config changed while it was being loaded")
    study = run_curvature_range_check(
        benchmark_config,
        check_config,
        progress=lambda message: print(message, flush=True),
    )
    return write_curvature_range_check_run(
        study,
        resolved,
        source_path,
        REPO_ROOT,
        generation_check_config_sha256=check_config_sha256,
        generation_benchmark_config_sha256=benchmark_config_sha256,
        generation_provenance=generation_provenance,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help=(
            "Path to one range-check config. The active ORCA contract has "
            "separate dual-port and dark-field configs under configs/."
        ),
    )
    arguments = parser.parse_args()
    for label, path in run(arguments.config).items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
