"""Camera binning and noise helpers extracted from the reference notebook."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import csr_matrix


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


def centered_camera_shape(
    input_shape: tuple[int, int],
    input_pixel_size_m: float,
    camera_pixel_size_m: float,
) -> tuple[int, int]:
    """Return the largest odd, centred camera array contained in the input grid.

    The numerical grids in this project place the origin at index ``n // 2``.
    An odd camera array is therefore used so that its central physical pixel is
    also centred on the origin. Only complete physical camera pixels are kept;
    any sub-pixel numerical margin at the field edge is excluded.
    """

    if len(input_shape) != 2 or any(size <= 0 for size in input_shape):
        raise ValueError("input_shape must contain two positive dimensions")
    if input_pixel_size_m <= 0 or camera_pixel_size_m <= 0:
        raise ValueError("pixel sizes must be positive")

    output: list[int] = []
    for size in input_shape:
        centres = (np.arange(size, dtype=float) - size // 2) * input_pixel_size_m
        left_edge = centres[0] - input_pixel_size_m / 2
        right_edge = centres[-1] + input_pixel_size_m / 2
        symmetric_width = 2 * min(-left_edge, right_edge)
        count = int(np.floor((symmetric_width + 1e-12 * camera_pixel_size_m) / camera_pixel_size_m))
        if count % 2 == 0:
            count -= 1
        if count <= 0:
            raise ValueError("camera pixel is larger than the centred numerical field")
        output.append(count)
    return output[0], output[1]


@lru_cache(maxsize=32)
def _centred_axis_weights(
    input_size: int,
    input_pixel_size_m: float,
    camera_pixel_size_m: float,
    output_size: int,
) -> csr_matrix:
    """Build one-dimensional overlap weights for physical-pixel integration."""

    input_centres = (
        np.arange(input_size, dtype=float) - input_size // 2
    ) * input_pixel_size_m
    input_left = input_centres - input_pixel_size_m / 2
    input_right = input_centres + input_pixel_size_m / 2

    output_centres = (
        np.arange(output_size, dtype=float) - (output_size - 1) / 2
    ) * camera_pixel_size_m
    output_left = output_centres - camera_pixel_size_m / 2
    output_right = output_centres + camera_pixel_size_m / 2

    if output_left[0] < input_left[0] - 1e-12 * input_pixel_size_m:
        raise ValueError("requested camera array extends beyond the numerical field")
    if output_right[-1] > input_right[-1] + 1e-12 * input_pixel_size_m:
        raise ValueError("requested camera array extends beyond the numerical field")

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for row, (left, right) in enumerate(zip(output_left, output_right, strict=True)):
        overlap = np.minimum(input_right, right) - np.maximum(input_left, left)
        indices = np.flatnonzero(overlap > 0)
        rows.extend([row] * len(indices))
        cols.extend(indices.tolist())
        values.extend((overlap[indices] / camera_pixel_size_m).tolist())

    weights = csr_matrix(
        (values, (rows, cols)),
        shape=(output_size, input_size),
        dtype=float,
    )
    row_sums = np.asarray(weights.sum(axis=1)).ravel()
    if not np.allclose(row_sums, 1.0, rtol=0.0, atol=2e-12):
        raise RuntimeError("physical camera-pixel overlap weights do not sum to unity")
    return weights


def resample_to_camera_pixels(
    image: ArrayLike,
    input_pixel_size_m: float,
    camera_pixel_size_m: float,
    output_shape: tuple[int, int] | None = None,
) -> NDArray:
    """Average a numerical image over centred physical camera pixels.

    Unlike :func:`bin_to_camera_pixels`, this operation does not require the
    physical pixel pitch to be an integer multiple of the numerical spacing.
    The input is treated as piecewise constant over each numerical cell and a
    separable overlap integral gives the mean intensity in every camera pixel.
    A constant image is preserved exactly and the output pixel area remains the
    physical area used by the photon-count calculation.
    """

    array = np.asarray(image)
    if array.ndim != 2:
        raise ValueError("camera resampling expects a 2D image")
    if input_pixel_size_m <= 0 or camera_pixel_size_m <= 0:
        raise ValueError("pixel sizes must be positive")
    if output_shape is None:
        output_shape = centered_camera_shape(
            array.shape,
            input_pixel_size_m,
            camera_pixel_size_m,
        )
    if len(output_shape) != 2 or any(size <= 0 for size in output_shape):
        raise ValueError("output_shape must contain two positive dimensions")

    row_weights = _centred_axis_weights(
        array.shape[0],
        float(input_pixel_size_m),
        float(camera_pixel_size_m),
        int(output_shape[0]),
    )
    column_weights = _centred_axis_weights(
        array.shape[1],
        float(input_pixel_size_m),
        float(camera_pixel_size_m),
        int(output_shape[1]),
    )
    row_resampled = row_weights @ array
    return np.asarray((column_weights @ row_resampled.T).T)


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
