from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import evaluate_faraday_operating_point


PARAMS = {
    "detuning_hz": 1.5e9,
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


REQUIRED_KEYS = {
    "faraday_signal_rad",
    "faraday_signal_scale",
    "scattered_photons_per_atom",
    "reabsorption_fraction",
    "destructiveness_metric",
    "estimated_per_frame_snr",
    "signal_per_scattered_photon",
    "information_per_scattered_photon",
    "signal_to_destruction",
}


def test_evaluate_faraday_operating_point_returns_required_finite_values() -> None:
    result = evaluate_faraday_operating_point(**PARAMS)

    assert REQUIRED_KEYS.issubset(result)
    for key in REQUIRED_KEYS:
        assert np.isfinite(result[key])


def test_evaluate_faraday_operating_point_non_negative_cost_metrics() -> None:
    result = evaluate_faraday_operating_point(**PARAMS)

    assert result["faraday_signal_scale"] >= 0
    assert result["scattered_photons_per_atom"] >= 0
    assert result["reabsorption_fraction"] >= 0
    assert result["destructiveness_metric"] >= 0
    assert result["signal_per_scattered_photon"] >= 0
    assert result["information_per_scattered_photon"] >= 0
    assert result["signal_to_destruction"] >= 0


def test_signal_per_scattered_photon_matches_simple_ratio() -> None:
    result = evaluate_faraday_operating_point(**PARAMS)

    expected = result["faraday_signal_scale"] / result["scattered_photons_per_atom"]
    assert result["signal_per_scattered_photon"] == pytest.approx(expected)
    assert result["information_per_scattered_photon"] == pytest.approx(expected)


def test_signal_per_scattered_photon_improves_when_power_is_reduced() -> None:
    lower_power = evaluate_faraday_operating_point(**{**PARAMS, "probe_power_mw": 1.0})
    higher_power = evaluate_faraday_operating_point(**{**PARAMS, "probe_power_mw": 2.0})

    assert lower_power["faraday_signal_scale"] == pytest.approx(higher_power["faraday_signal_scale"])
    assert lower_power["scattered_photons_per_atom"] < higher_power["scattered_photons_per_atom"]
    assert lower_power["signal_per_scattered_photon"] > higher_power["signal_per_scattered_photon"]


def test_evaluate_faraday_operating_point_is_deterministic() -> None:
    first = evaluate_faraday_operating_point(**PARAMS)
    second = evaluate_faraday_operating_point(**PARAMS)

    assert first == second
