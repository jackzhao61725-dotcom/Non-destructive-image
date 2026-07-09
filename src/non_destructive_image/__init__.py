"""Shared helpers for the non-destructive imaging notebook refactor.

The reference notebook remains authoritative. These helpers are mechanical
extractions of repeated notebook patterns and are intentionally small so they can
be regression-checked before being wired back into the notebook exports.
"""

from .atomic_model import ThomasFermiState, build_thomas_fermi_state, recoil_quantities
from .analysis import evaluate_faraday_operating_point, sweep_faraday_detuning
from .camera import (
    add_camera_noise,
    bin_to_camera_pixels,
    normalize_camera_counts,
    simulate_camera_image,
    simulate_noisy_camera_image,
)
from .fourier import propagate_scattered_field
from .imaging import simulate_dgi_image, simulate_faraday_image, simulate_fourier_image, simulate_pci_image
from .light_atom import (
    dimensionless_detuning,
    faraday_rotation_angle,
    intensity_at_atoms,
    reabsorption_fraction,
    residual_optical_depth,
    scalar_phase_shift,
    scattered_photons_per_atom,
)
from .multishot import accumulate_snr, simulate_multishot_sequence
from .profiles import thomas_fermi_profile_2d

__all__ = [
    "ThomasFermiState",
    "accumulate_snr",
    "add_camera_noise",
    "bin_to_camera_pixels",
    "build_thomas_fermi_state",
    "evaluate_faraday_operating_point",
    "normalize_camera_counts",
    "simulate_camera_image",
    "simulate_noisy_camera_image",
    "dimensionless_detuning",
    "faraday_rotation_angle",
    "intensity_at_atoms",
    "reabsorption_fraction",
    "propagate_scattered_field",
    "recoil_quantities",
    "residual_optical_depth",
    "scalar_phase_shift",
    "scattered_photons_per_atom",
    "simulate_fourier_image",
    "simulate_dgi_image",
    "simulate_faraday_image",
    "simulate_pci_image",
    "simulate_multishot_sequence",
    "sweep_faraday_detuning",
    "thomas_fermi_profile_2d",
]
