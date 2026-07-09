from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import sweep_faraday_intensity


PARAMS = {
    "detuning_hz": 1.5e9,
    "column_density_peak": 2.0e14,
    "resonant_cross_section": 7.68e-14,
    "gamma_rad_per_s": 2 * np.pi * 29.5e6,
    "pulse_duration_s": 40e-6,
    "saturation_intensity": 560.0,
    "probe_diameter_m": 24e-3,
    "kappa_f": 1.0,
    "column_densities_for_reabsorption": np.array([2.0e14, 1.0e14, 1.5e14]),
    "photons_per_camera_pixel": 500.0,
}

PROBE_POWERS_MW = np.array([1.0, 2.0, 4.0])

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


def test_faraday_intensity_sweep_returns_finite_arrays_with_expected_length() -> None:
    result = sweep_faraday_intensity(PROBE_POWERS_MW, **PARAMS)

    np.testing.assert_allclose(result["probe_power_mw"], PROBE_POWERS_MW)
    for key in ARRAY_KEYS:
        assert result[key].shape == PROBE_POWERS_MW.shape
        assert np.all(np.isfinite(result[key]))


def test_faraday_intensity_sweep_best_index_matches_objective_maximum() -> None:
    result = sweep_faraday_intensity(PROBE_POWERS_MW, **PARAMS)

    metric = result["signal_per_scattered_photon"]
    expected_best = int(np.argmax(metric))
    assert result["objective_key"] == "signal_per_scattered_photon"
    assert result["best_index"] == expected_best
    assert result["best_probe_power_mw"] == pytest.approx(PROBE_POWERS_MW[expected_best])
    assert result["best_objective_value"] == pytest.approx(metric[expected_best])


def test_faraday_intensity_sweep_is_deterministic() -> None:
    first = sweep_faraday_intensity(PROBE_POWERS_MW, **PARAMS)
    second = sweep_faraday_intensity(PROBE_POWERS_MW, **PARAMS)

    for key in {"probe_power_mw", *ARRAY_KEYS}:
        np.testing.assert_allclose(first[key], second[key])
    assert first["objective_key"] == second["objective_key"]
    assert first["best_index"] == second["best_index"]
    assert first["best_probe_power_mw"] == pytest.approx(second["best_probe_power_mw"])
    assert first["best_objective_value"] == pytest.approx(second["best_objective_value"])


def test_faraday_intensity_sweep_accepts_single_probe_power() -> None:
    result = sweep_faraday_intensity([2.0], **PARAMS)

    assert result["probe_power_mw"].shape == (1,)
    assert result["best_index"] == 0
    assert result["best_probe_power_mw"] == pytest.approx(2.0)


def test_faraday_intensity_sweep_rejects_empty_input() -> None:
    with pytest.raises(ValueError, match="non-empty 1D array"):
        sweep_faraday_intensity([], **PARAMS)


@pytest.mark.parametrize("bad_probe_powers", ([0.0, 1.0], [-1.0, 1.0]))
def test_faraday_intensity_sweep_rejects_non_positive_probe_power(
    bad_probe_powers: list[float],
) -> None:
    with pytest.raises(ValueError, match="positive values"):
        sweep_faraday_intensity(bad_probe_powers, **PARAMS)
