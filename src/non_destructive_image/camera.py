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


def simulate_camera_image(
    image: ArrayLike,
    bin_size: int = 15,
    photons_per_pixel: float | None = None,
    *,
    return_intermediates: bool = False,
) -> NDArray[np.floating] | dict[str, NDArray[np.floating]]:
    """Run the deterministic notebook camera pipeline without stochastic noise.

    The notebook camera path first bins the high-resolution image to camera
    pixels. When a photon scale is supplied, this helper also applies the
    deterministic count conversion and normalisation used around the stochastic
    camera model, but intentionally does not add Poisson or read noise.
    """

    binned_image = bin_to_camera_pixels(image, bin_size)
    if photons_per_pixel is None:
        if return_intermediates:
            return {"binned_image": binned_image, "camera_image": binned_image}
        return binned_image

    deterministic_counts = binned_image * photons_per_pixel
    camera_image = normalize_camera_counts(deterministic_counts, photons_per_pixel)
    if return_intermediates:
        return {
            "binned_image": binned_image,
            "deterministic_counts": deterministic_counts,
            "camera_image": camera_image,
        }
    return camera_image


def simulate_noisy_camera_image(
    image: ArrayLike,
    photons_per_pixel: float,
    rng: np.random.Generator,
    read_noise_electrons: float,
    bin_size: int = 15,
    *,
    input_is_binned: bool = False,
    normalize: bool = True,
    return_intermediates: bool = False,
) -> NDArray[np.floating] | dict[str, NDArray[np.floating]]:
    """Run the notebook stochastic camera recipe with explicit RNG handling.

    This helper is a thin orchestration layer around the existing camera
    helpers. It optionally bins a high-resolution image, applies
    ``add_camera_noise(...)`` using the caller-provided random generator, and
    optionally normalises the noisy electron counts back to image units.
    """

    binned_image = np.asarray(image) if input_is_binned else bin_to_camera_pixels(image, bin_size)
    noisy_counts = add_camera_noise(
        binned_image,
        photons_per_pixel,
        rng,
        read_noise_electrons,
    )
    noisy_image = normalize_camera_counts(noisy_counts, photons_per_pixel) if normalize else noisy_counts

    if return_intermediates:
        return {
            "binned_image": binned_image,
            "noisy_counts": noisy_counts,
            "noisy_image": noisy_image,
        }
    return noisy_image
