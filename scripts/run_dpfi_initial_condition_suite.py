"""Run the frozen DPFI inverse on the declared initial-condition suite."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies.initial_condition_suite import (
    run_and_write_initial_condition_suite,
    validate_initial_condition_suite,
)
from non_destructive_image.reconstruction.studies.io import load_json


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    REPO_ROOT / "configs" / "dpfi_initial_condition_suite_v1_orca_fusion_m10.json"
)


def validate_only(config_path: Path) -> dict[str, object]:
    """Validate source lineage, frozen candidate and all analytic truth maps."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = load_json(resolved)
    context = validate_initial_condition_suite(config, REPO_ROOT)
    return {
        "source_run_id": context.source.run_id,
        "condition_count": len(context.truths),
        "condition_ids": tuple(truth.morphology.name for truth in context.truths),
        "fluence_mw_us": context.fluences_mw_us,
        "realizations_per_condition": context.realizations_per_condition,
        "expected_trial_count": (
            len(context.truths)
            * len(context.fluences_mw_us)
            * context.realizations_per_condition
        ),
        "integration_grid_shape": context.integration_support.shape,
        "supported_cell_count": int(context.integration_support.support_mask.sum()),
        "candidate": context.candidate.label,
        "parameter_count": context.candidate.model.parameter_count,
    }


def run(config_path: Path) -> dict[str, Path]:
    """Run and seal the configured DPFI initial-condition sweep."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = load_json(resolved)
    return run_and_write_initial_condition_suite(
        config,
        resolved,
        REPO_ROOT,
        progress=lambda message: print(message, flush=True),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="validate the suite contract and truth maps without fitting or writing",
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
