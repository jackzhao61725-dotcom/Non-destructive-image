import json
import csv
from pathlib import Path

import numpy as np

from scripts.generate_accumulated_snr_invariance_plot import build_accumulated_snr_data


ROOT = Path(__file__).resolve().parents[2]


def _configs():
    with (ROOT / "configs" / "dissertation_plots_v1.json").open(encoding="utf-8") as handle:
        plot = json.load(handle)
    with (ROOT / plot["notebook_defaults_config"]).open(encoding="utf-8") as handle:
        notebook = json.load(handle)
    return notebook, plot


def test_accumulated_snr_scaling_and_mode_independent_budget() -> None:
    data, params, checks = build_accumulated_snr_data(*_configs())

    assert params["destruction_budget_fraction"] == 0.30
    assert checks["passed"] is True
    assert checks["n_max_mode_independent"] is True
    assert 1.9 <= checks["n_max_log_slope"] <= 2.1
    assert -2.2 <= checks["quadratic_shot_plus_read_snr_high_detuning_log_slope"] <= -0.9
    assert -1.1 <= checks["quadratic_shot_noise_only_snr_shot_log_slope"] <= -0.9
    assert checks["pci_shot_plus_read_total_relative_range"] <= 0.01
    assert checks["pci_shot_noise_total_relative_range"] <= 0.01
    assert checks["quadratic_shot_noise_total_relative_range"] <= 0.03
    assert checks["pci_shot_plus_read_le_shot_noise_limit"] is True
    assert checks["dgi_shot_plus_read_le_shot_noise_limit"] is True
    assert checks["pci_and_dgi_shot_noise_curves_coincide"] is False
    assert 1.9 <= checks["pci_to_dgi_shot_noise_total_ratio_median"] <= 2.1
    assert checks["pci_to_dgi_shot_noise_ratio_relative_range"] <= 0.04
    assert checks["notebook_linear_phase_regime_satisfied"] is True
    assert checks["full_fourier_camera_pipeline_used"] is False
    assert np.all(data["pci_shot_plus_read_snr_total"] <= data["pci_shot_only_snr_total"])
    assert np.all(data["dgi_shot_plus_read_snr_total"] <= data["dgi_shot_only_snr_total"])
    assert data["dgi_shot_plus_read_snr_total"][-1] < data["dgi_shot_plus_read_snr_total"][0]
    assert np.all(np.isfinite(data["n_max"]))
    assert np.count_nonzero(np.isclose(data["detuning_ghz"], 1.5, rtol=0.0, atol=1e-12)) == 1


def test_plot_output_contract_exists_after_generation() -> None:
    output = ROOT / "results" / "dissertation_plots_v1" / "accumulated_snr_invariance"
    assert (output / "accumulated_snr_invariance.svg").exists()
    assert (output / "figure_3_2.png").exists()
    assert (output / "figure_3_2.pdf").exists()
    assert (output / "accumulated_snr_invariance_data.csv").exists()
    assert (output / "figure_3_2_parameter_register.csv").exists()
    assert (output / "accumulated_snr_invariance_summary.json").exists()
    assert (output / "metadata.json").exists()


def test_figure_3_2_parameter_register_has_required_provenance_columns() -> None:
    path = ROOT / "results" / "dissertation_plots_v1" / "accumulated_snr_invariance" / "figure_3_2_parameter_register.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 5
    assert set(rows[0]) == {
        "量", "值", "Δ/2π (GHz)", "P (mW)", "τ (µs)", "imaging axis", "归一化",
        "N_max模型", "QE/read", "repo 路径",
    }
    assert all(float(row["Δ/2π (GHz)"]) == 1.5 for row in rows)
    assert all(float(row["P (mW)"]) == 1.0 for row in rows)
    assert all(float(row["τ (µs)"]) == 90.0 for row in rows)
    assert all(row["imaging axis"] == "x (0)" for row in rows)
