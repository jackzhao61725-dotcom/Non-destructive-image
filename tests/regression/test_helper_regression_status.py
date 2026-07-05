from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import (
    build_thomas_fermi_state,
    dimensionless_detuning,
    faraday_rotation_angle,
    propagate_scattered_field,
    reabsorption_fraction,
    residual_optical_depth,
    scalar_phase_shift,
    thomas_fermi_profile_2d,
)


def test_core_helper_regression_against_original_formulas() -> None:
    hbar = 1.054571817e-34
    k_b = 1.380649e-23
    amu = 1.66053907e-27
    a0 = 5.29177211e-11
    atomic_mass = 166 * amu
    scattering_length = 72 * a0
    trap_hz = np.array([293.0, 14.0, 233.0])
    state = build_thomas_fermi_state(2.5e4, scattering_length, trap_hz, atomic_mass, hbar, k_b)

    expected_omega = 2 * np.pi * trap_hz
    expected_omega_bar = expected_omega.prod() ** (1 / 3)
    expected_a_ho = np.sqrt(hbar / (atomic_mass * expected_omega_bar))
    expected_mu = 0.5 * (15 * 2.5e4 * scattering_length / expected_a_ho) ** (2 / 5) * hbar * expected_omega_bar
    np.testing.assert_allclose(state.radii, np.sqrt(2 * expected_mu / (atomic_mass * expected_omega**2)))

    detuning_hz = 1.5e9
    gamma = 2 * np.pi * 29.5e6
    sigma0 = 3 * (401.0e-9) ** 2 / (2 * np.pi)
    delta = 2 * detuning_hz * 2 * np.pi / gamma
    assert dimensionless_detuning(detuning_hz, gamma) == pytest.approx(delta)
    assert scalar_phase_shift(detuning_hz, state.column_density[0], sigma0, gamma) == pytest.approx(
        sigma0 * state.column_density[0] * delta / (2 * (1 + delta**2))
    )
    assert residual_optical_depth(detuning_hz, state.column_density[0], sigma0, gamma) == pytest.approx(
        sigma0 * state.column_density[0] / (1 + delta**2)
    )
    assert faraday_rotation_angle(detuning_hz, state.column_density[0], sigma0, gamma, 1.0) == pytest.approx(
        sigma0 * state.column_density[0] * delta / (2 * (1 + delta**2))
    )
    expected_od = sigma0 * state.column_density / (1 + delta**2)
    assert reabsorption_fraction(detuning_hz, state.column_density, sigma0, gamma) == pytest.approx(
        float(np.mean(1 - np.exp(-expected_od)))
    )

    axis = np.linspace(-2e-6, 2e-6, 8)
    grid_a, grid_b = np.meshgrid(axis, axis)
    profile = thomas_fermi_profile_2d(grid_a, grid_b, state.radii[0], state.radii[2])
    np.testing.assert_allclose(profile, np.maximum(0, 1 - grid_a**2 / state.radii[0] ** 2 - grid_b**2 / state.radii[2] ** 2) ** 1.5)

    pupil = np.ones_like(profile)
    scattered = np.exp(1j * 0.1 * profile) - 1
    np.testing.assert_allclose(propagate_scattered_field(scattered, pupil), np.fft.ifft2(np.fft.fft2(scattered) * pupil))
