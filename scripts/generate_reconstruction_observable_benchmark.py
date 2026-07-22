"""Replay and seal the frozen held-out observable reconstruction benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies.io import load_json
from non_destructive_image.reconstruction.studies.observable_benchmark import (
    build_observable_integration_support,
    run_and_write_observable_benchmark,
    validate_source_benchmark_artifacts,
)
from non_destructive_image.reconstruction.studies.morphology import (
    build_morphology_study_context,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    REPO_ROOT
    / "configs"
    / "reconstruction_observables_v1_orca_fusion_m10.json"
)


def run(config_path: Path) -> dict[str, Path]:
    """Validate the source benchmark, replay only held-out fits, and seal outputs."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = load_json(resolved)
    return run_and_write_observable_benchmark(
        config,
        resolved,
        REPO_ROOT,
        progress=lambda message: print(message, flush=True),
    )


def validate_only(config_path: Path) -> dict[str, object]:
    """Validate source lineage, held-out axes and support without running a fit."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = load_json(resolved)
    source = validate_source_benchmark_artifacts(config, REPO_ROOT)
    context = build_morphology_study_context(source.config)
    support = build_observable_integration_support(config, source, context)
    return {
        "source_run_id": source.run_id,
        "held_out_trial_count": len(source.held_out_rows),
        "integration_grid_shape": support.shape,
        "supported_cell_count": int(support.support_mask.sum()),
        "integration_physical_area_m2": support.physical_area_m2,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="validate the frozen source and support without replaying any fit",
    )
    arguments = parser.parse_args()
    outputs = (
        validate_only(arguments.config)
        if arguments.validate_only
        else run(arguments.config)
    )
    for label, value in outputs.items():
        print(f"{label}: {value}")


if __name__ == "__main__":
    main()
