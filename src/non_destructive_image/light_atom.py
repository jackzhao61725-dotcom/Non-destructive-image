"""Light-atom interaction formulas extracted from the reference notebook.

All functions take physical constants and experimental parameters explicitly so
parameter values remain controlled by the notebook/reference configuration.
"""

from __future__ import annotations

import numpy as np


def dimensionless_detuning(detuning_hz: float, gamma_rad_per_s: float) -> float:
    """Return the notebook detuning parameter ``delta = 2 * Delta / Gamma``.

    The notebook stores ``Gamma`` in rad/s and ``Delta`` in Hz, so the exact
    implementation is ``2 * detuning_hz * 2*pi / gamma_rad_per_s``.
    """

    return 2 * detuning_hz * 2 * np.pi / gamma_rad_per_s


def scalar_phase_shift(
    detuning_hz: float,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
) -> float:
    """Return the scalar dispersive phase shift used by the notebook."""

    detuning = dimensionless_detuning(detuning_hz, gamma_rad_per_s)
    return resonant_cross_section * column_density_peak * detuning / (2 * (1 + detuning**2))


def residual_optical_depth(
    detuning_hz: float,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
) -> float:
    """Return the residual on-resonance-scaled optical depth at detuning."""

    detuning = dimensionless_detuning(detuning_hz, gamma_rad_per_s)
    return resonant_cross_section * column_density_peak / (1 + detuning**2)


def intensity_at_atoms(
    probe_power_mw: float,
    probe_diameter_m: float,
    use_peak_intensity: bool = True,
) -> float:
    """Convert probe power to intensity at the cloud using notebook semantics.

    The reference notebook uses twice the area-averaged Gaussian intensity when
    ``use_peak_intensity`` is true because the atoms sit at the beam centre.
    """

    area_averaged = (probe_power_mw * 1e-3) / (np.pi * (probe_diameter_m / 2) ** 2)
    return 2 * area_averaged if use_peak_intensity else area_averaged


def scattered_photons_per_atom(
    detuning_hz: float,
    probe_power_mw: float,
    pulse_duration_s: float,
    saturation_intensity: float,
    gamma_rad_per_s: float,
    probe_diameter_m: float,
    use_peak_intensity: bool = True,
) -> float:
    """Return scattered photons per atom per shot using the notebook formula."""

    saturation_parameter = intensity_at_atoms(
        probe_power_mw,
        probe_diameter_m,
        use_peak_intensity=use_peak_intensity,
    ) / saturation_intensity
    detuning = dimensionless_detuning(detuning_hz, gamma_rad_per_s)
    return (gamma_rad_per_s / 2) * saturation_parameter / (1 + saturation_parameter + detuning**2) * pulse_duration_s


def faraday_rotation_angle(
    detuning_hz: float,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
    kappa_f: float,
) -> float:
    """Return the signed peak Faraday rotation ``theta_F = kappa_F * phi``.

    ``kappa_f`` is supplied by the atomic-response model. For the fully
    spin-polarised axial 166Er reference state used in the dissertation it is
    ``-45/91``; apparatus-level polarimetric calibration is handled separately.
    """

    return kappa_f * scalar_phase_shift(
        detuning_hz,
        column_density_peak,
        resonant_cross_section,
        gamma_rad_per_s,
    )


def reabsorption_fraction(
    detuning_hz: float,
    column_densities: np.ndarray,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
) -> float:
    """Return the notebook's angle-averaged Rayleigh reabsorption fraction.

    This is the direct extraction of ``mean(1 - exp(-OD))`` with
    ``OD = sigma0 * n_col / (1 + delta**2)`` along the three principal axes.
    """

    detuning = dimensionless_detuning(detuning_hz, gamma_rad_per_s)
    optical_depth = (
        resonant_cross_section * np.asarray(column_densities) / (1 + detuning**2)
    )
    return float(np.mean(1 - np.exp(-optical_depth)))
