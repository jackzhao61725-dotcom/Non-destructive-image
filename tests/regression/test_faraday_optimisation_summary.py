from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import summarise_faraday_sweep, sweep_faraday_detuning


PARAMS = {
    "column_density_peak": 2.0e14,
    "resonant_cross_section": 7.68e-14,
    "gamma_rad_per_s": 2 * np.pi * 29.5e6,
    "probe_power_mw": 2.0,
    "pulse_duration_s": 40e-6,
    "saturation_intensity": 560.0,
    "probe_diameter_m": 24e-3,
    "kappa_f": 1.0,
    "column_densities_for_reabsorption": np.array([2.0e14, 1.0e14, 1.5e14]),
    "photons_per_camera_pixel": 500.0,
}

DETUNINGS = np.array([1.0e9, 1.5e9, 2.0e9])


def test_faraday_sweep_summary_best_index_matches_metric_maximum() -> None:
    sweep = sweep_faraday_detuning(DETUNINGS, **PARAMS)
    summary = summarise_faraday_sweep(sweep, "signal_per_scattered_photon")

    metric = sweep["signal_per_scattered_photon"]
    expected_best = int(np.argmax(metric))
    assert summary["metric_key"] == "signal_per_scattered_photon"
    assert summary["parameter_key"] == "detuning_hz"
    assert summary["best_index"] == expected_best
    assert summary["best_parameter_value"] == pytest.approx(DETUNINGS[expected_best])
    assert summary["best_metric_value"] == pytest.approx(metric[expected_best])
    assert summary["num_points"] == len(DETUNINGS)
    assert summary["metric_min"] == pytest.approx(np.min(metric))
    assert summary["metric_max"] == pytest.approx(np.max(metric))


def test_faraday_sweep_summary_handles_single_point_sweep() -> None:
    sweep = sweep_faraday_detuning([1.5e9], **PARAMS)
    summary = summarise_faraday_sweep(sweep)

    assert summary["best_index"] == 0
    assert summary["best_parameter_value"] == pytest.approx(1.5e9)
    assert summary["num_points"] == 1


def test_faraday_sweep_summary_rejects_missing_metric() -> None:
    sweep = sweep_faraday_detuning(DETUNINGS, **PARAMS)

    with pytest.raises(ValueError, match="metric_key is not available"):
        summarise_faraday_sweep(sweep, "missing_metric")


def test_faraday_sweep_summary_rejects_empty_sweep() -> None:
    empty_sweep = {
        "detuning_hz": np.array([]),
        "signal_per_scattered_photon": np.array([]),
    }

    with pytest.raises(ValueError, match="non-empty 1D array"):
        summarise_faraday_sweep(empty_sweep)


def test_faraday_sweep_summary_is_deterministic() -> None:
    sweep = sweep_faraday_detuning(DETUNINGS, **PARAMS)

    first = summarise_faraday_sweep(sweep)
    second = summarise_faraday_sweep(sweep)

    assert first == second
