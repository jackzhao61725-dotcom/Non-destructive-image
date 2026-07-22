"""Run and serialize the frozen morphology reconstruction benchmark."""

from __future__ import annotations

import argparse
from pathlib import Path

from non_destructive_image.reconstruction.studies import (
    file_sha256,
    load_json,
    run_morphology_benchmark_study,
    write_morphology_benchmark_run,
)
from non_destructive_image.reconstruction.studies.provenance import (
    capture_reconstruction_provenance,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = (
    REPO_ROOT
    / "configs"
    / "reconstruction_morphology_benchmark_v4_orca_fusion_m10.json"
)


def run(config_path: Path) -> dict[str, Path]:
    """Execute one declared study and write its immutable artifacts."""

    resolved = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    generation_provenance = capture_reconstruction_provenance(
        REPO_ROOT,
        entry_points=(Path("scripts/generate_reconstruction_morphology_benchmark.py"),),
    )
    generation_config_sha256 = file_sha256(resolved)
    config = load_json(resolved)
    if file_sha256(resolved) != generation_config_sha256:
        raise RuntimeError("benchmark config changed while it was being loaded")
    study = run_morphology_benchmark_study(config, progress=lambda message: print(message, flush=True))
    return write_morphology_benchmark_run(
        study,
        resolved,
        REPO_ROOT,
        generation_config_sha256=generation_config_sha256,
        generation_provenance=generation_provenance,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    arguments = parser.parse_args()
    for label, path in run(arguments.config).items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
