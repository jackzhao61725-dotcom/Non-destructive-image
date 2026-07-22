"""Render reconstruction benchmark figures from frozen artifacts only."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies import (
    generate_morphology_benchmark_figures,
    load_json,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    REPO_ROOT
    / "configs"
    / "reconstruction_morphology_benchmark_v4_orca_fusion_m10.json"
)


def run(config_path: Path) -> dict[str, Path]:
    """Resolve the output location and render only the recorded artifacts."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = load_json(resolved)
    configured = Path(str(config["output_directory"]))
    output_directory = configured if configured.is_absolute() else REPO_ROOT / configured
    return generate_morphology_benchmark_figures(output_directory)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    arguments = parser.parse_args()
    for label, path in run(arguments.config).items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
