"""Small deterministic analysis helpers built on the migrated core."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from .light_atom import faraday_rotation_angle, reabsorption_fraction, scattered_photons_per_atom


def evaluate_faraday_operating_point(
    detuning_hz: float,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
    probe_power_mw: float,
    pulse_duration_s: float,
    saturation_intensity: float,
    probe_diameter_m: float,
    *,
    kappa_f: float = 1.0,
    column_densities_for_reabsorption: ArrayLike | None = None,
    use_peak_intensity: bool = True,
    photons_per_camera_pixel: float | None = None,
) -> dict[str, float]:
    """Evaluate one deterministic Faraday operating point.

    This helper is intentionally small: it combines existing notebook-equivalent
    light-atom helpers into scalar quantities useful for later optimisation. It
    does not scan parameters, add stochastic noise, plot results, or alter the
    current phenomenological ``kappa_f`` convention.
    """

    faraday_signal = faraday_rotation_angle(
        detuning_hz,
        column_density_peak,
        resonant_cross_section,
        gamma_rad_per_s,
        kappa_f,
    )
    scattered_photons = scattered_photons_per_atom(
        detuning_hz,
        probe_power_mw,
        pulse_duration_s,
        saturation_intensity,
        gamma_rad_per_s,
        probe_diameter_m,
        use_peak_intensity=use_peak_intensity,
    )
    reabsorption = (
        reabsorption_fraction(
            detuning_hz,
            np.asarray(column_densities_for_reabsorption),
            resonant_cross_section,
            gamma_rad_per_s,
        )
        if column_densities_for_reabsorption is not None
        else np.nan
    )

    destructiveness = scattered_photons * (1 + (0.0 if np.isnan(reabsorption) else reabsorption))
    signal_scale = abs(faraday_signal)
    signal_per_scattered_photon = signal_scale / scattered_photons if scattered_photons > 0 else np.inf
    signal_to_destruction = signal_scale / destructiveness if destructiveness > 0 else np.inf
    estimated_per_frame_snr = (
        2 * signal_scale * np.sqrt(photons_per_camera_pixel)
        if photons_per_camera_pixel is not None
        else np.nan
    )

    return {
        "faraday_signal_rad": float(faraday_signal),
        "faraday_signal_scale": float(signal_scale),
        "scattered_photons_per_atom": float(scattered_photons),
        "reabsorption_fraction": float(reabsorption),
        "destructiveness_metric": float(destructiveness),
        "estimated_per_frame_snr": float(estimated_per_frame_snr),
        "signal_per_scattered_photon": float(signal_per_scattered_photon),
        "information_per_scattered_photon": float(signal_per_scattered_photon),
        "signal_to_destruction": float(signal_to_destruction),
    }


def sweep_faraday_detuning(
    detuning_hz_values: ArrayLike,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
    probe_power_mw: float,
    pulse_duration_s: float,
    saturation_intensity: float,
    probe_diameter_m: float,
    *,
    kappa_f: float = 1.0,
    column_densities_for_reabsorption: ArrayLike | None = None,
    use_peak_intensity: bool = True,
    photons_per_camera_pixel: float | None = None,
    objective_key: str = "signal_per_scattered_photon",
) -> dict[str, np.ndarray | float | int | str]:
    """Evaluate a deterministic one-dimensional Faraday detuning sweep."""

    detunings = np.asarray(detuning_hz_values, dtype=float)
    if detunings.ndim != 1 or detunings.size == 0:
        raise ValueError("detuning_hz_values must be a non-empty 1D array")

    operating_points = [
        evaluate_faraday_operating_point(
            float(detuning_hz),
            column_density_peak,
            resonant_cross_section,
            gamma_rad_per_s,
            probe_power_mw,
            pulse_duration_s,
            saturation_intensity,
            probe_diameter_m,
            kappa_f=kappa_f,
            column_densities_for_reabsorption=column_densities_for_reabsorption,
            use_peak_intensity=use_peak_intensity,
            photons_per_camera_pixel=photons_per_camera_pixel,
        )
        for detuning_hz in detunings
    ]
    if objective_key not in operating_points[0]:
        raise ValueError(f"objective_key is not available: {objective_key}")

    result_arrays = {
        key: np.asarray([point[key] for point in operating_points], dtype=float)
        for key in operating_points[0]
    }
    objective_values = result_arrays[objective_key]
    if np.all(np.isnan(objective_values)):
        raise ValueError(f"objective_key contains only NaN values: {objective_key}")

    best_index = int(np.nanargmax(objective_values))
    return {
        "detuning_hz": detunings,
        "objective_key": objective_key,
        "best_index": best_index,
        "best_detuning_hz": float(detunings[best_index]),
        "best_objective_value": float(objective_values[best_index]),
        **result_arrays,
    }
