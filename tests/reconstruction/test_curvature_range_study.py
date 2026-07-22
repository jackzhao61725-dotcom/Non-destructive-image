from __future__ import annotations

from pathlib import Path

import pytest

from non_destructive_image.reconstruction.studies.curvature_range import (
    build_curvature_range_check_context,
)
from non_destructive_image.reconstruction.studies.io import load_json


REPO_ROOT = Path(__file__).resolve().parents[2]
ORCA_BENCHMARK_CONFIG = (
    REPO_ROOT
    / "configs"
    / "reconstruction_morphology_benchmark_v4_orca_fusion_m10.json"
)
ORCA_CHECK_CONFIGS = (
    REPO_ROOT
    / "configs"
    / "reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port.json",
    REPO_ROOT
    / "configs"
    / "reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field.json",
)


@pytest.mark.parametrize("check_path", ORCA_CHECK_CONFIGS)
def test_active_range_checks_change_only_declared_basis_and_weights(
    check_path: Path,
) -> None:
    benchmark = load_json(ORCA_BENCHMARK_CONFIG)
    context = build_curvature_range_check_context(benchmark, load_json(check_path))

    assert [candidate.label for candidate in context.candidates] == [
        "resolution_matched_17x5__curvature_30_um2",
        "resolution_matched_17x5__curvature_100_um2",
        "resolution_matched_17x5__curvature_300_um2",
        "resolution_matched_17x5__curvature_1000_um2",
    ]
    assert all(candidate.model.parameter_count == 85 for candidate in context.candidates)
    assert context.config["ensemble"] == benchmark["ensemble"]
    assert context.config["fit"] == benchmark["fit"]
    assert context.config["detector"] == benchmark["detector"]
    assert context.config["physics"] == benchmark["physics"]
    assert [case.name for case in context.morphology_split.calibration] == benchmark[
        "synthetic_morphologies"
    ]["calibration_names"]
    assert context.grid.sampling_mode == "physical_pixel"
    assert context.grid.camera_shape == (153, 153)
    assert context.reference_count_scale_e == pytest.approx(220.58087277528466)
    assert context.read_noise_e == pytest.approx(1.4)
