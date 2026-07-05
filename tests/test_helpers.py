from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import (
    add_camera_noise,
    bin_to_camera_pixels,
    normalize_camera_counts,
    propagate_scattered_field,
    simulate_fourier_image,
    thomas_fermi_profile_2d,
)


def test_thomas_fermi_profile_matches_notebook_expression() -> None:
    axis = np.linspace(-2.0, 2.0, 7)
    grid_a, grid_b = np.meshgrid(axis, axis)
    radius_a = 1.5
    radius_b = 0.75

    expected = np.maximum(0, 1 - grid_a**2 / radius_a**2 - grid_b**2 / radius_b**2) ** 1.5

    np.testing.assert_allclose(thomas_fermi_profile_2d(grid_a, grid_b, radius_a, radius_b), expected)


def test_propagate_scattered_field_matches_notebook_fft_pattern() -> None:
    rng = np.random.default_rng(123)
    scattered = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    pupil = (rng.random((8, 8)) > 0.25).astype(float)

    expected = np.fft.ifft2(np.fft.fft2(scattered) * pupil)

    np.testing.assert_allclose(propagate_scattered_field(scattered, pupil), expected)



def test_simulate_fourier_image_matches_notebook_pci_dgi_core() -> None:
    axis = np.linspace(-1.0, 1.0, 8)
    grid_a, grid_b = np.meshgrid(axis, axis)
    profile = np.maximum(0, 1 - grid_a**2 / 0.9**2 - grid_b**2 / 0.7**2) ** 1.5
    phase_peak = 0.2
    object_field = np.exp(1j * phase_peak * profile)
    pupil = ((grid_a**2 + grid_b**2) < 0.85).astype(float)

    scattered = np.exp(1j * phase_peak * profile) - 1
    propagated = np.fft.ifft2(np.fft.fft2(scattered) * pupil)

    pci_reference = 0.95 * np.exp(1j * np.pi / 2)
    dgi_reference = 10 ** (-4.0 / 2)

    np.testing.assert_allclose(
        simulate_fourier_image(object_field, pupil, pci_reference),
        np.abs(pci_reference + propagated) ** 2,
    )
    np.testing.assert_allclose(
        simulate_fourier_image(object_field, pupil, dgi_reference),
        np.abs(dgi_reference + propagated) ** 2,
    )
    np.testing.assert_allclose(
        simulate_fourier_image(object_field, pupil, pci_reference, return_intensity=False),
        pci_reference + propagated,
    )

def test_bin_to_camera_pixels_matches_notebook_reshape_mean() -> None:
    image = np.arange(31 * 46, dtype=float).reshape(31, 46)
    bin_size = 15
    rows = (image.shape[0] // bin_size) * bin_size
    cols = (image.shape[1] // bin_size) * bin_size
    expected = image[:rows, :cols].reshape(rows // bin_size, bin_size, cols // bin_size, bin_size).mean(axis=(1, 3))

    np.testing.assert_allclose(bin_to_camera_pixels(image, bin_size), expected)


def test_camera_noise_matches_notebook_recipe_for_seeded_rng() -> None:
    binned = np.array([[1.0, -0.5], [0.25, 0.0]])
    photons_per_pixel = 10.0
    read_noise = 0.5

    expected_rng = np.random.default_rng(7)
    expected = expected_rng.poisson(np.clip(binned, 0, None) * photons_per_pixel) + expected_rng.normal(
        0,
        read_noise,
        binned.shape,
    )

    actual_rng = np.random.default_rng(7)
    np.testing.assert_allclose(add_camera_noise(binned, photons_per_pixel, actual_rng, read_noise), expected)


def test_normalize_camera_counts_divides_by_photons_per_pixel() -> None:
    counts = np.array([[2.0, 4.0]])

    np.testing.assert_allclose(normalize_camera_counts(counts, 2.0), np.array([[1.0, 2.0]]))

from non_destructive_image import (
    dimensionless_detuning,
    faraday_rotation_angle,
    intensity_at_atoms,
    reabsorption_fraction,
    residual_optical_depth,
    scalar_phase_shift,
    scattered_photons_per_atom,
)


def test_light_atom_helpers_match_notebook_formulas() -> None:
    detuning_hz = 1.5e9
    gamma = 2 * np.pi * 29.5e6
    column_density = 1.2e14
    sigma0 = 7.68e-14
    detuning = 2 * detuning_hz * 2 * np.pi / gamma

    assert dimensionless_detuning(detuning_hz, gamma) == pytest.approx(detuning)
    assert scalar_phase_shift(detuning_hz, column_density, sigma0, gamma) == pytest.approx(
        sigma0 * column_density * detuning / (2 * (1 + detuning**2))
    )
    assert residual_optical_depth(detuning_hz, column_density, sigma0, gamma) == pytest.approx(
        sigma0 * column_density / (1 + detuning**2)
    )

    kappa_f = 1.0
    assert faraday_rotation_angle(
        detuning_hz, column_density, sigma0, gamma, kappa_f
    ) == pytest.approx(
        kappa_f * sigma0 * column_density * detuning / (2 * (1 + detuning**2))
    )

    column_densities = np.array(
        [column_density, 0.5 * column_density, 0.25 * column_density]
    )
    expected_od = sigma0 * column_densities / (1 + detuning**2)
    assert reabsorption_fraction(detuning_hz, column_densities, sigma0, gamma) == pytest.approx(
        float(np.mean(1 - np.exp(-expected_od)))
    )


def test_intensity_and_scattering_helpers_match_notebook_formulas() -> None:
    probe_power_mw = 2.0
    probe_diameter_m = 24e-3
    gamma = 2 * np.pi * 29.5e6
    saturation_intensity = 600.0
    detuning_hz = 1.5e9
    pulse_duration_s = 40e-6

    area_averaged = (probe_power_mw * 1e-3) / (np.pi * (probe_diameter_m / 2) ** 2)
    peak_intensity = 2 * area_averaged
    saturation_parameter = peak_intensity / saturation_intensity
    detuning = 2 * detuning_hz * 2 * np.pi / gamma

    assert intensity_at_atoms(probe_power_mw, probe_diameter_m, use_peak_intensity=True) == pytest.approx(
        peak_intensity
    )
    assert intensity_at_atoms(probe_power_mw, probe_diameter_m, use_peak_intensity=False) == pytest.approx(
        area_averaged
    )
    assert scattered_photons_per_atom(
        detuning_hz,
        probe_power_mw,
        pulse_duration_s,
        saturation_intensity,
        gamma,
        probe_diameter_m,
    ) == pytest.approx((gamma / 2) * saturation_parameter / (1 + saturation_parameter + detuning**2) * pulse_duration_s)

from non_destructive_image import build_thomas_fermi_state, recoil_quantities


def test_build_thomas_fermi_state_matches_notebook_algebra() -> None:
    hbar = 1.054571817e-34
    k_b = 1.380649e-23
    amu = 1.66053907e-27
    a0 = 5.29177211e-11
    atom_number = 2.5e4
    atomic_mass = 166 * amu
    scattering_length = 72 * a0
    trap_hz = np.array([293.0, 14.0, 233.0])

    omega = 2 * np.pi * trap_hz
    omega_bar = omega.prod() ** (1 / 3)
    a_ho = np.sqrt(hbar / (atomic_mass * omega_bar))
    mu = 0.5 * (15 * atom_number * scattering_length / a_ho) ** (2 / 5) * hbar * omega_bar
    n_peak = mu * atomic_mass / (4 * np.pi * hbar**2 * scattering_length)
    radii = np.sqrt(2 * mu / (atomic_mass * omega**2))
    column_density = (4 / 3) * n_peak * radii
    atom_number_check = (8 * np.pi / 15) * n_peak * radii.prod()

    state = build_thomas_fermi_state(atom_number, scattering_length, trap_hz, atomic_mass, hbar, k_b)

    np.testing.assert_allclose(state.trap_angular_frequencies, omega)
    assert state.geometric_mean_frequency == pytest.approx(omega_bar)
    assert state.harmonic_oscillator_length == pytest.approx(a_ho)
    assert state.chemical_potential == pytest.approx(mu)
    assert state.chemical_potential_temperature == pytest.approx(mu / k_b)
    assert state.peak_density == pytest.approx(n_peak)
    np.testing.assert_allclose(state.radii, radii)
    np.testing.assert_allclose(state.column_density, column_density)
    assert state.atom_number_check == pytest.approx(atom_number_check)


def test_recoil_quantities_match_notebook_definitions() -> None:
    hbar = 1.054571817e-34
    k_b = 1.380649e-23
    atomic_mass = 166 * 1.66053907e-27
    wavevector = 2 * np.pi / 401e-9

    expected_energy = (hbar * wavevector) ** 2 / (2 * atomic_mass)
    expected_temperature = expected_energy / k_b
    expected_velocity = hbar * wavevector / atomic_mass

    recoil_energy, recoil_temperature, recoil_velocity = recoil_quantities(hbar, wavevector, atomic_mass, k_b)

    assert recoil_energy == pytest.approx(expected_energy)
    assert recoil_temperature == pytest.approx(expected_temperature)
    assert recoil_velocity == pytest.approx(expected_velocity)
