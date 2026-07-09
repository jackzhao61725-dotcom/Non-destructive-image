from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import sweep_faraday_exposure_time


PARAMS = {
    "detuning_hz": 1.5e9,
    "column_density_peak": 2.0e14,
    "resonant_cross_section": 7.68e-14,
    "gamma_rad_per_s": 2 * np.pi * 29.5e6,
    "probe_power_mw": 2.0,
    "saturation_intensity": 560.0,
    "probe_diameter_m": 24e-3,
    "kappa_f": 1.0,
    "column_densities_for_reabsorption": np.array([2.0e14, 1.0e14, 1.5e14]),
    "photons_per_camera_pixel": 500.0,
}

PULSE_DURATIONS_S = np.array([20e-6, 40e-6, 80e-6])

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


def test_faraday_exposure_sweep_returns_finite_arrays_with_expected_length() -> None:
    result = sweep_faraday_exposure_time(PULSE_DURATIONS_S, **PARAMS)

    np.testing.assert_allclose(result["pulse_duration_s"], PULSE_DURATIONS_S)
    for key in ARRAY_KEYS:
        assert result[key].shape == PULSE_DURATIONS_S.shape
        assert np.all(np.isfinite(result[key]))


def test_faraday_exposure_sweep_best_index_matches_objective_maximum() -> None:
    result = sweep_faraday_exposure_time(PULSE_DURATIONS_S, **PARAMS)

    metric = result["signal_per_scattered_photon"]
    expected_best = int(np.argmax(metric))
    assert result["objective_key"] == "signal_per_scattered_photon"
    assert result["best_index"] == expected_best
    assert result["best_pulse_duration_s"] == pytest.approx(PULSE_DURATIONS_S[expected_best])
    assert result["best_objective_value"] == pytest.approx(metric[expected_best])


def test_faraday_exposure_sweep_is_deterministic() -> None:
    first = sweep_faraday_exposure_time(PULSE_DURATIONS_S, **PARAMS)
    second = sweep_faraday_exposure_time(PULSE_DURATIONS_S, **PARAMS)

    for key in {"pulse_duration_s", *ARRAY_KEYS}:
        np.testing.assert_allclose(first[key], second[key])
    assert first["objective_key"] == second["objective_key"]
    assert first["best_index"] == second["best_index"]
    assert first["best_pulse_duration_s"] == pytest.approx(second["best_pulse_duration_s"])
    assert first["best_objective_value"] == pytest.approx(second["best_objective_value"])


def test_faraday_exposure_sweep_accepts_single_exposure_time() -> None:
    result = sweep_faraday_exposure_time([40e-6], **PARAMS)

    assert result["pulse_duration_s"].shape == (1,)
    assert result["best_index"] == 0
    assert result["best_pulse_duration_s"] == pytest.approx(40e-6)


def test_faraday_exposure_sweep_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="non-empty 1D array"):
        sweep_faraday_exposure_time([], **PARAMS)


@pytest.mark.parametrize("bad_pulse_durations", ([0.0, 40e-6], [-20e-6, 40e-6]))
def test_faraday_exposure_sweep_rejects_non_positive_exposure_time(
    bad_pulse_durations: list[float],
) -> None:
    with pytest.raises(ValueError, match="positive values"):
        sweep_faraday_exposure_time(bad_pulse_durations, **PARAMS)
