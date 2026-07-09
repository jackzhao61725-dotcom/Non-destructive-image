from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_condensate_stage import build_condensate_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_condensate_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_condensate_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_profile"].shape == (1024, 1024)
    assert stage["notebook_column_density_m2"].shape == (1024, 1024)
    assert report["profile_comparison"]["max_absolute_difference"] == 0.0
    assert report["column_density_map_comparison"]["max_absolute_difference"] == 0.0
    assert report["vector_comparisons"]["radii_m"]["max_absolute_difference"] == 0.0
    assert report["vector_comparisons"]["column_density_m2"]["max_absolute_difference"] == 0.0


def test_notebook_condensate_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_condensate_stage(config)
    state = stage["state"]

    np.testing.assert_allclose(
        state["radii"],
        np.array([1.1857996475343457e-06, 2.4817092623397375e-05, 1.4911557799466233e-06]),
        rtol=1e-12,
        atol=0.0,
    )
    np.testing.assert_allclose(
        state["column_density"],
        np.array([5.3759624525784675e14, 1.1251121418610648e16, 6.760330466117985e14]),
        rtol=1e-12,
        atol=0.0,
    )
    assert stage["notebook_profile"][512, 512] == 1.0
    assert stage["notebook_column_density_m2"][512, 512] == pytest.approx(state["column_density"][0])
    assert np.unravel_index(np.argmax(stage["notebook_column_density_m2"]), (1024, 1024)) == (512, 512)
