from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import sweep_faraday_detuning


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

ARRAY_KEYS = {
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


def test_faraday_detuning_sweep_returns_finite_arrays_with_expected_length() -> None:
    result = sweep_faraday_detuning(DETUNINGS, **PARAMS)

    np.testing.assert_allclose(result["detuning_hz"], DETUNINGS)
    for key in ARRAY_KEYS:
        assert result[key].shape == DETUNINGS.shape
        assert np.all(np.isfinite(result[key]))


def test_faraday_detuning_sweep_best_index_matches_objective_maximum() -> None:
    result = sweep_faraday_detuning(DETUNINGS, **PARAMS)

    metric = result["signal_per_scattered_photon"]
    expected_best = int(np.argmax(metric))
    assert result["objective_key"] == "signal_per_scattered_photon"
    assert result["best_index"] == expected_best
    assert result["best_detuning_hz"] == pytest.approx(DETUNINGS[expected_best])
    assert result["best_objective_value"] == pytest.approx(metric[expected_best])


def test_faraday_detuning_sweep_is_deterministic() -> None:
    first = sweep_faraday_detuning(DETUNINGS, **PARAMS)
    second = sweep_faraday_detuning(DETUNINGS, **PARAMS)

    for key in {"detuning_hz", *ARRAY_KEYS}:
        np.testing.assert_allclose(first[key], second[key])
    assert first["objective_key"] == second["objective_key"]
    assert first["best_index"] == second["best_index"]
    assert first["best_detuning_hz"] == pytest.approx(second["best_detuning_hz"])
    assert first["best_objective_value"] == pytest.approx(second["best_objective_value"])


def test_faraday_detuning_sweep_accepts_single_detuning() -> None:
    result = sweep_faraday_detuning([1.5e9], **PARAMS)

    assert result["detuning_hz"].shape == (1,)
    assert result["best_index"] == 0
    assert result["best_detuning_hz"] == pytest.approx(1.5e9)


def test_faraday_detuning_sweep_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="non-empty 1D array"):
        sweep_faraday_detuning([], **PARAMS)
