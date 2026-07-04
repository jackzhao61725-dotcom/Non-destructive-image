"""Camera binning and noise helpers extracted from the reference notebook."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def bin_to_camera_pixels(image: ArrayLike, bin_size: int = 15) -> NDArray[np.floating]:
    """Bin a high-resolution image to camera pixels using notebook semantics.

    The notebook truncates each axis to an integer multiple of 15 grid cells and
    averages over each 15-by-15 block. The default ``bin_size`` preserves that
    convention while making the repeated calculation explicit.
    """

    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError("camera binning expects a 2D image")
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")

    rows = (image.shape[0] // bin_size) * bin_size
    cols = (image.shape[1] // bin_size) * bin_size
    trimmed = image[:rows, :cols]
    return trimmed.reshape(rows // bin_size, bin_size, cols // bin_size, bin_size).mean(axis=(1, 3))


def add_camera_noise(
    binned_image: ArrayLike,
    photons_per_pixel: float,
    rng: np.random.Generator,
    read_noise_electrons: float,
) -> NDArray[np.floating]:
    """Add Poisson photon noise and Gaussian read noise in notebook units.

    This preserves the notebook recipe:
    ``poisson(clip(image, 0, None) * photons_per_pixel) + normal(0, read_noise)``.
    The returned values are electron counts, matching the intermediate notebook
    variable before division by ``photons_per_pixel``.
    """

    binned_image = np.asarray(binned_image)
    return rng.poisson(np.clip(binned_image, 0, None) * photons_per_pixel) + rng.normal(
        0,
        read_noise_electrons,
        binned_image.shape,
    )


def normalize_camera_counts(counts: ArrayLike, photons_per_pixel: float) -> NDArray[np.floating]:
    """Convert noisy camera counts back to the notebook's normalised image units."""

    return np.asarray(counts) / photons_per_pixel
