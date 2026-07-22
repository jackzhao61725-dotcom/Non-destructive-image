"""Shared coherent Fourier-imaging helpers for PCI/DGI-style paths."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .fourier import propagate_scattered_field


def simulate_fourier_image(
    object_field: ArrayLike,
    pupil: ArrayLike,
    reference_field: complex | float,
    *,
    return_intensity: bool = True,
) -> NDArray[np.floating] | NDArray[np.complexfloating]:
    """Propagate the notebook's scattered field and recombine with a reference.

    The PCI/DGI notebook path propagates only the scattered component
    ``object_field - 1`` through ``ifft2(fft2(scattered) * pupil)`` and then
    adds a mode-specific unscattered reference field before taking ``abs(E)**2``.
    This helper preserves the notebook FFT convention, normalisation, phase sign,
    and mask multiplication order.
    """

    scattered_field = np.asarray(object_field) - 1
    image_field = reference_field + propagate_scattered_field(scattered_field, pupil)
    if return_intensity:
        return np.abs(image_field) ** 2
    return image_field


def simulate_pci_image(
    phase_map: ArrayLike,
    pupil: ArrayLike,
    phase_plate_transmittance: float = 0.95,
    phase_plate_phase: float = np.pi / 2,
    *,
    return_intermediates: bool = False,
) -> NDArray[np.floating] | dict[str, NDArray[np.floating] | NDArray[np.complexfloating] | complex]:
    """Return the notebook-equivalent scalar PCI image intensity.

    This is the PCI-specific orchestration above ``simulate_fourier_image``:
    ``object_field = exp(1j * phase_map)`` and reference field
    ``t_p * exp(1j * theta)``. It preserves the notebook phase-plate,
    scattered-field, FFT/pupil, and intensity conventions.
    """

    object_field = np.exp(1j * np.asarray(phase_map))
    pci_reference_field = phase_plate_transmittance * np.exp(1j * phase_plate_phase)
    image_field = simulate_fourier_image(
        object_field,
        pupil,
        pci_reference_field,
        return_intensity=False,
    )
    pci_image_intensity = np.abs(image_field) ** 2

    if return_intermediates:
        return {
            "object_field": object_field,
            "scattered_field": object_field - 1,
            "propagated_scattered_field": image_field - pci_reference_field,
            "pci_reference_field": pci_reference_field,
            "pci_image_intensity": pci_image_intensity,
        }
    return pci_image_intensity


def simulate_dgi_image(
    phase_map: ArrayLike,
    pupil: ArrayLike,
    stop_optical_depth: float = 4.0,
    *,
    return_intermediates: bool = False,
) -> NDArray[np.floating] | dict[str, NDArray[np.floating] | NDArray[np.complexfloating] | float]:
    """Return the notebook-equivalent scalar DGI image intensity.

    This is the DGI-specific orchestration above ``simulate_fourier_image``:
    ``object_field = exp(1j * phase_map)`` and attenuated carrier reference
    ``10**(-OD/2)``. It preserves the notebook scattered-field, FFT/pupil,
    dark-ground attenuation, and intensity conventions.
    """

    object_field = np.exp(1j * np.asarray(phase_map))
    dgi_reference_field = 10 ** (-stop_optical_depth / 2)
    image_field = simulate_fourier_image(
        object_field,
        pupil,
        dgi_reference_field,
        return_intensity=False,
    )
    dgi_image_intensity = np.abs(image_field) ** 2

    if return_intermediates:
        return {
            "object_field": object_field,
            "scattered_field": object_field - 1,
            "propagated_scattered_field": image_field - dgi_reference_field,
            "dgi_reference_field": dgi_reference_field,
            "dgi_image_intensity": dgi_image_intensity,
        }
    return dgi_image_intensity


def simulate_faraday_image(
    theta_f_map: ArrayLike,
    pupil: ArrayLike,
    *,
    return_intermediates: bool = False,
) -> dict[str, NDArray[np.floating] | NDArray[np.complexfloating]]:
    """Return Faraday dark-field and signed dual-port outputs.

    The caller supplies the signed rotation map ``theta_F``. The helper applies
    opposite circular phase shifts, the common FFT/pupil propagation, and the
    circular-to-linear recombination. Historical ``u``/``v`` output keys are
    retained for regression compatibility; explicit analyser ``H``/``V``
    aliases expose the dissertation convention
    ``S = (I_H - I_V) / (I_H + I_V)``.
    """

    theta_f_map_array = np.asarray(theta_f_map)
    sigma_plus_object_field = np.exp(1j * theta_f_map_array)
    sigma_minus_object_field = np.exp(-1j * theta_f_map_array)

    sigma_plus_scattered_field = sigma_plus_object_field - 1
    sigma_minus_scattered_field = sigma_minus_object_field - 1
    sigma_plus_propagated_scattered_field = propagate_scattered_field(
        sigma_plus_scattered_field,
        pupil,
    )
    sigma_minus_propagated_scattered_field = propagate_scattered_field(
        sigma_minus_scattered_field,
        pupil,
    )
    sigma_plus_field = 1 + sigma_plus_propagated_scattered_field
    sigma_minus_field = 1 + sigma_minus_propagated_scattered_field

    output_ex_field = (sigma_plus_field + sigma_minus_field) / 2
    output_ey_field = 1j * (sigma_plus_field - sigma_minus_field) / 2
    dark_field_intensity = np.abs(output_ey_field) ** 2
    dual_port_u_intensity = np.abs(output_ex_field + output_ey_field) ** 2 / 2
    dual_port_v_intensity = np.abs(output_ex_field - output_ey_field) ** 2 / 2
    dual_port_signal = (dual_port_v_intensity - dual_port_u_intensity) / (
        dual_port_v_intensity + dual_port_u_intensity
    )

    outputs = {
        "dark_field_intensity": dark_field_intensity,
        "dual_port_u_intensity": dual_port_u_intensity,
        "dual_port_v_intensity": dual_port_v_intensity,
        "analyser_h_intensity": dual_port_v_intensity,
        "analyser_v_intensity": dual_port_u_intensity,
        "dual_port_signal": dual_port_signal,
    }
    if return_intermediates:
        return {
            "theta_f_map_rad": theta_f_map_array,
            "sigma_plus_object_field": sigma_plus_object_field,
            "sigma_minus_object_field": sigma_minus_object_field,
            "sigma_plus_scattered_field": sigma_plus_scattered_field,
            "sigma_minus_scattered_field": sigma_minus_scattered_field,
            "sigma_plus_propagated_scattered_field": sigma_plus_propagated_scattered_field,
            "sigma_minus_propagated_scattered_field": sigma_minus_propagated_scattered_field,
            "sigma_plus_field": sigma_plus_field,
            "sigma_minus_field": sigma_minus_field,
            "output_ex_field": output_ex_field,
            "output_ey_field": output_ey_field,
            "parallel_field": output_ex_field,
            "perpendicular_field": -output_ey_field,
            **outputs,
        }
    return outputs
