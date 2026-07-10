from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.generate_condensate_three_view import build_three_view_stage, three_view_summary


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_condensate_three_view_contains_all_projection_axes() -> None:
    config = _load_config()
    stage = build_three_view_stage(config)

    assert [view["integrated_axis_label"] for view in stage["views"]] == ["x", "y", "z"]
    assert [view["display_plane_labels"] for view in stage["views"]] == [["y", "z"], ["x", "z"], ["x", "y"]]
    for view in stage["views"]:
        assert view["column_density_m2"].shape == (1024, 1024)
        assert view["absolute_column_density"] is True
        assert view["normalised"] is False
        assert np.all(np.isfinite(view["column_density_m2"]))
        assert np.min(view["column_density_m2"]) >= 0.0


def test_condensate_three_view_peaks_match_tf_column_densities() -> None:
    config = _load_config()
    stage = build_three_view_stage(config)
    expected = stage["condensate_stage"]["state"]["column_density"]

    for axis, view in enumerate(stage["views"]):
        image = view["column_density_m2"]
        assert image[512, 512] == pytest.approx(expected[axis], rel=1e-12)
        assert np.max(image) == pytest.approx(expected[axis], rel=1e-12)
        assert np.unravel_index(np.argmax(image), image.shape) == (512, 512)


def test_condensate_three_view_summary_records_model_extension_boundary() -> None:
    config = _load_config()
    stage = build_three_view_stage(config)
    summary = three_view_summary(config, stage)

    assert summary["not_exact_notebook_figure"] is True
    assert summary["no_experimental_calibration_applied"] is True
    assert summary["not_final_calibrated_prediction"] is True
    assert summary["quantity"] == "absolute Thomas-Fermi column density"
    assert summary["normalised"] is False
    assert len(summary["views"]) == 3
