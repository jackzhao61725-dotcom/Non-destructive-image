from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import bin_to_camera_pixels, normalize_camera_counts, simulate_camera_image


def test_simulate_camera_image_matches_existing_binning_helper() -> None:
    image = np.arange(31 * 46, dtype=float).reshape(31, 46)
    bin_size = 15

    result = simulate_camera_image(image, bin_size=bin_size)
    expected = bin_to_camera_pixels(image, bin_size=bin_size)

    assert result.shape == (2, 3)
    assert np.isfinite(result).all()
    np.testing.assert_allclose(result, expected)


def test_simulate_camera_image_applies_deterministic_normalisation() -> None:
    image = np.linspace(0.25, 1.25, 30 * 45).reshape(30, 45)
    bin_size = 15
    photons_per_pixel = 100.0

    result = simulate_camera_image(
        image,
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        return_intermediates=True,
    )
    expected_binned = bin_to_camera_pixels(image, bin_size=bin_size)
    expected_counts = expected_binned * photons_per_pixel
    expected_camera = normalize_camera_counts(expected_counts, photons_per_pixel)

    assert result["camera_image"].shape == (2, 3)
    assert np.isfinite(result["camera_image"]).all()
    np.testing.assert_allclose(result["binned_image"], expected_binned)
    np.testing.assert_allclose(result["deterministic_counts"], expected_counts)
    np.testing.assert_allclose(result["camera_image"], expected_camera)
    np.testing.assert_allclose(result["camera_image"], expected_binned)
