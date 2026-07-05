"""Fourier-optics helpers extracted from the reference notebook."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def propagate_scattered_field(scattered_field: ArrayLike, pupil: ArrayLike) -> NDArray[np.complexfloating]:
    """Apply the notebook FFT -> pupil -> inverse FFT propagation pattern.

    The reference notebook repeatedly propagates only the scattered part of the
    field, e.g. ``exp(1j * phase_map) - 1``, using
    ``ifft2(fft2(scattered_field) * pupil)``. This helper preserves that exact
    convention and deliberately does not apply shifts, padding, or alternative
    Fourier normalisations.
    """

    return np.fft.ifft2(np.fft.fft2(scattered_field) * np.asarray(pupil))
