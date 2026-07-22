from __future__ import annotations

from scripts.generate_reconstruction_morphology_benchmark import (
    DEFAULT_CONFIG as BENCHMARK_CONFIG,
)
from scripts.plot_reconstruction_morphology_benchmark import (
    DEFAULT_CONFIG as BENCHMARK_PLOT_CONFIG,
)
from scripts.run_reconstruction_credibility_study import (
    DEFAULT_CONFIG as CREDIBILITY_CONFIG,
)
from scripts.run_reconstruction_curvature_range_check import ACTIVE_CONFIGS


def test_default_reconstruction_entry_points_use_the_active_orca_contract() -> None:
    assert BENCHMARK_CONFIG.name == (
        "reconstruction_morphology_benchmark_v4_orca_fusion_m10.json"
    )
    assert BENCHMARK_PLOT_CONFIG == BENCHMARK_CONFIG
    assert CREDIBILITY_CONFIG.name == (
        "reconstruction_credibility_v2_orca_fusion_m10.json"
    )


def test_active_curvature_checks_require_both_faraday_readouts() -> None:
    assert {path.name for path in ACTIVE_CONFIGS} == {
        "reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port.json",
        "reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field.json",
    }
