"""Run the reconstruction credibility study without using truth error as evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies.credibility import (
    run_and_write_credibility_study,
)
from non_destructive_image.reconstruction.studies.io import load_json


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    REPO_ROOT / "configs" / "reconstruction_credibility_v2_orca_fusion_m10.json"
)


def run(config_path: Path) -> dict[str, Path]:
    """Load one declared study and serialize its credibility artifacts."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    study_config = load_json(resolved)
    source_relative = Path(str(study_config["source_benchmark_config"]))
    source_path = (
        source_relative
        if source_relative.is_absolute()
        else REPO_ROOT / source_relative
    )
    source_config = load_json(source_path)
    return run_and_write_credibility_study(
        study_config,
        resolved,
        source_config,
        source_path,
        REPO_ROOT,
        progress=lambda message: print(message, flush=True),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    arguments = parser.parse_args()
    for label, path in run(arguments.config).items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
