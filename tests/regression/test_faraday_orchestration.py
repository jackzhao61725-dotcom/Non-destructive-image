from __future__ import annotations

from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import simulate_faraday_image


BASELINE_PATH = Path("regression/baseline/imaging/faraday_imaging_baseline_v1.npz")


def test_simulate_faraday_image_matches_faraday_baseline() -> None:
    with np.load(BASELINE_PATH) as baseline:
        result = simulate_faraday_image(
            baseline["theta_f_map_rad"],
            baseline["pupil"],
            return_intermediates=True,
        )

        for key in [
            "theta_f_map_rad",
            "sigma_plus_object_field",
            "sigma_minus_object_field",
            "sigma_plus_field",
            "sigma_minus_field",
            "output_ex_field",
            "output_ey_field",
            "dark_field_intensity",
            "dual_port_u_intensity",
            "dual_port_v_intensity",
            "dual_port_signal",
        ]:
            np.testing.assert_allclose(result[key], baseline[key], rtol=1e-10, atol=1e-12)


def test_reversing_signed_rotation_swaps_analyser_ports() -> None:
    theta_f = np.linspace(-0.13, 0.11, 99).reshape(9, 11)
    pupil = np.ones_like(theta_f)

    forward = simulate_faraday_image(theta_f, pupil, return_intermediates=True)
    reversed_rotation = simulate_faraday_image(-theta_f, pupil, return_intermediates=True)

    np.testing.assert_allclose(
        reversed_rotation["sigma_plus_field"],
        forward["sigma_minus_field"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reversed_rotation["sigma_minus_field"],
        forward["sigma_plus_field"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reversed_rotation["analyser_h_intensity"],
        forward["analyser_v_intensity"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reversed_rotation["analyser_v_intensity"],
        forward["analyser_h_intensity"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reversed_rotation["dark_field_intensity"],
        forward["dark_field_intensity"],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        reversed_rotation["dual_port_signal"],
        -forward["dual_port_signal"],
        rtol=1e-12,
        atol=1e-12,
    )


def test_negative_erbium_kappa_gives_negative_signal_in_h_minus_v_convention() -> None:
    phase = np.full((9, 9), 0.203)
    theta_f = (-45 / 91) * phase
    result = simulate_faraday_image(theta_f, np.ones_like(theta_f))

    np.testing.assert_allclose(
        result["analyser_h_intensity"],
        result["dual_port_v_intensity"],
        rtol=0,
        atol=0,
    )
    np.testing.assert_allclose(
        result["analyser_v_intensity"],
        result["dual_port_u_intensity"],
        rtol=0,
        atol=0,
    )
    np.testing.assert_allclose(
        result["dual_port_signal"],
        np.sin(2 * theta_f),
        rtol=1e-12,
        atol=1e-12,
    )
    assert np.all(result["dual_port_signal"] < 0)
    assert np.all(result["analyser_h_intensity"] < result["analyser_v_intensity"])
