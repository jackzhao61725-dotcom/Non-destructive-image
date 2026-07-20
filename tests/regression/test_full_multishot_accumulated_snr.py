import csv
import json
from pathlib import Path

import numpy as np

from scripts.generate_full_multishot_accumulated_snr import (
    matched_dual_port_snr,
    matched_filter_snr,
)


ROOT = Path(__file__).resolve().parents[2]


def test_matched_filter_snr_orders_shot_and_read_noise() -> None:
    image = np.array([[1.0, 1.1], [0.9, 1.2]])
    roi = np.ones_like(image, dtype=bool)
    result = matched_filter_snr(image, 1.0, 1000.0, 7.0, roi)
    assert result["snr_shot_plus_read_noise"] <= result["snr_shot_noise_only"]
    assert result["signal_l2_e"] > 0
    assert np.isfinite(result["shot_plus_read_effective_e"])


def test_dual_port_matched_snr_uses_sum_variance_and_two_readouts() -> None:
    port_u = np.full((2, 2), 0.45)
    port_v = np.full((2, 2), 0.55)
    roi = np.ones_like(port_u, dtype=bool)
    result = matched_dual_port_snr(port_u, port_v, 1000.0, 7.0, roi)

    signal_per_pixel = 100.0
    expected_shot = np.sqrt(4 * signal_per_pixel**2 / 1000.0)
    expected_shot_read = np.sqrt(4 * signal_per_pixel**2 / (1000.0 + 2 * 7.0**2))
    np.testing.assert_allclose(result["snr_shot_noise_only"], expected_shot)
    np.testing.assert_allclose(result["snr_shot_plus_read_noise"], expected_shot_read)


def test_generated_full_multishot_summary_contract() -> None:
    output = ROOT / "results" / "dissertation_plots_v1" / "full_multishot_accumulated_snr"
    required = [
        "full_multishot_accumulated_snr.svg",
        "full_multishot_accumulated_snr_data.csv",
        "framewise_snr_sequences.csv",
        "model_comparison.csv",
        "full_multishot_accumulated_snr_summary.json",
        "metadata.json",
        "faraday_canonical_reference_at_1p5GHz.csv",
    ]
    assert all((output / name).exists() for name in required)
    with (output / "full_multishot_accumulated_snr_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)
    assert summary["checks"]["passed"] is True
    assert summary["checks"]["pci_shot_plus_read_le_shot_all"] is True
    assert summary["checks"]["dgi_shot_plus_read_le_shot_all"] is True
    assert summary["checks"]["faraday_dark_field_shot_plus_read_le_shot_all"] is True
    assert summary["checks"]["faraday_dual_port_shot_plus_read_le_shot_all"] is True
    assert summary["checks"]["faraday_kappa_f_is_v1_placeholder"] is True
    assert set(summary["canonical_faraday_reference_1p5GHz"]) == {
        "Faraday dark-field",
        "Faraday dual-port",
    }
    for values in summary["canonical_faraday_reference_1p5GHz"].values():
        assert values["shot_noise_only"] > 0
        assert 0 < values["shot_plus_read_noise"] <= values["shot_noise_only"]
    expected = {
        "Faraday dark-field": {
            "shot_noise_only": 55.77833731293387,
            "shot_plus_read_noise": 34.70877083139104,
        },
        "Faraday dual-port": {
            "shot_noise_only": 107.48197734510573,
            "shot_plus_read_noise": 106.555526352218,
        },
    }
    for mode, values in expected.items():
        for noise_model, value in values.items():
            np.testing.assert_allclose(
                summary["canonical_faraday_reference_1p5GHz"][mode][noise_model],
                value,
                rtol=1e-12,
                atol=1e-12,
            )
    assert "not calibrated absolute Faraday predictions" in summary["faraday_interpretation"]
    assert summary["verdict"].startswith("A.")

    with (output / "faraday_canonical_reference_at_1p5GHz.csv").open(
        encoding="utf-8",
        newline="",
    ) as handle:
        ledger = list(csv.DictReader(handle))
    assert len(ledger) == 4
    assert {int(row["roi_pixel_count"]) for row in ledger} == {228}
    assert {int(row["n_frames_full"]) for row in ledger} == {10}
    assert {float(row["initial_atom_number"]) for row in ledger} == {25000.0}
    assert {float(row["kappa_F"]) for row in ledger} == {1.0}
    assert {float(row["probe_power_mW"]) for row in ledger} == {1.0}
    assert all(np.isclose(float(row["exposure_time_us"]), 90.0) for row in ledger)
    assert {row["imaging_axis"] for row in ledger} == {"x"}
    assert all("uncalibrated" in row["calibration_status"] for row in ledger)
