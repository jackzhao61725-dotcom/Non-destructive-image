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
