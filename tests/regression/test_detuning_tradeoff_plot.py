from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")


OUTPUT_DIR = Path("results/dissertation_plots_v2_orca_fusion/detuning_tradeoff")
DATA_PATH = OUTPUT_DIR / "detuning_tradeoff_data.csv"
SVG_PATH = OUTPUT_DIR / "detuning_tradeoff.svg"
METADATA_PATH = OUTPUT_DIR / "metadata.json"


def _load_column(name: str) -> np.ndarray:
    with DATA_PATH.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return np.asarray([float(row[name]) for row in rows], dtype=float)


def _log_slope(x: np.ndarray, y: np.ndarray) -> float:
    count = max(4, int(np.ceil(x.size * 0.35)))
    tail_x = x[-count:]
    tail_y = y[-count:]
    return float(np.polyfit(np.log(tail_x), np.log(tail_y), 1)[0])


def test_detuning_tradeoff_outputs_exist() -> None:
    assert DATA_PATH.exists()
    assert SVG_PATH.exists()
    assert METADATA_PATH.exists()
    assert SVG_PATH.stat().st_size > 1000


def test_detuning_tradeoff_detuning_is_monotonic() -> None:
    detuning = _load_column("detuning_hz")

    assert detuning.ndim == 1
    assert detuning.size > 10
    assert np.all(np.diff(detuning) > 0)


def test_detuning_tradeoff_far_detuned_scaling() -> None:
    detuning = _load_column("detuning_hz")
    phase = _load_column("abs_scalar_phase_rad")
    scattering = _load_column("scattered_photons_per_atom")

    assert _log_slope(detuning, phase) == pytest.approx(-1.0, abs=0.04)
    assert _log_slope(detuning, scattering) == pytest.approx(-2.0, abs=0.04)


def test_detuning_tradeoff_metadata_records_provenance() -> None:
    with METADATA_PATH.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    plot_config = metadata["config_files_used"]["plot_config"].replace("\\", "/")
    assert plot_config == "configs/dissertation_plots_v2_orca_fusion.json"
    assert "dimensionless delta" in metadata["detuning_convention"]
    assert metadata["normalisation"]["absolute_values_saved"] is True
    assert metadata["calibration_status"] == "No experimental RAI / absorption calibration has been applied."
    assert metadata["prediction_status"].endswith("not a final operating-point prediction.")
