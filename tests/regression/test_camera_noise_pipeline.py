from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import (
    add_camera_noise,
    bin_to_camera_pixels,
    normalize_camera_counts,
    simulate_noisy_camera_image,
)


def test_simulate_noisy_camera_image_is_reproducible_with_same_seed() -> None:
    image = np.linspace(0.1, 1.2, 30 * 45).reshape(30, 45)

    first = simulate_noisy_camera_image(
        image,
        photons_per_pixel=250.0,
        rng=np.random.default_rng(7),
        read_noise_electrons=3.0,
    )
    second = simulate_noisy_camera_image(
        image,
        photons_per_pixel=250.0,
        rng=np.random.default_rng(7),
        read_noise_electrons=3.0,
    )

    np.testing.assert_allclose(first, second)


def test_simulate_noisy_camera_image_changes_with_different_seed() -> None:
    image = np.linspace(0.1, 1.2, 30 * 45).reshape(30, 45)

    first = simulate_noisy_camera_image(
        image,
        photons_per_pixel=250.0,
        rng=np.random.default_rng(7),
        read_noise_electrons=3.0,
    )
    second = simulate_noisy_camera_image(
        image,
        photons_per_pixel=250.0,
        rng=np.random.default_rng(8),
        read_noise_electrons=3.0,
    )

    assert not np.array_equal(first, second)


def test_simulate_noisy_camera_image_shape_and_finite_values() -> None:
    image = np.linspace(0.0, 1.0, 31 * 46).reshape(31, 46)
    result = simulate_noisy_camera_image(
        image,
        photons_per_pixel=100.0,
        rng=np.random.default_rng(11),
        read_noise_electrons=1.0,
        bin_size=15,
        return_intermediates=True,
    )

    assert result["binned_image"].shape == (2, 3)
    assert result["noisy_counts"].shape == (2, 3)
    assert result["noisy_image"].shape == (2, 3)
    assert np.isfinite(result["noisy_counts"]).all()
    assert np.isfinite(result["noisy_image"]).all()


def test_simulate_noisy_camera_image_matches_direct_helper_composition() -> None:
    image = np.linspace(0.1, 1.3, 30 * 45).reshape(30, 45)
    photons_per_pixel = 125.0
    read_noise_electrons = 2.5
    seed = 123

    expected_binned = bin_to_camera_pixels(image, bin_size=15)
    expected_counts = add_camera_noise(
        expected_binned,
        photons_per_pixel,
        np.random.default_rng(seed),
        read_noise_electrons,
    )
    expected_image = normalize_camera_counts(expected_counts, photons_per_pixel)

    result = simulate_noisy_camera_image(
        image,
        photons_per_pixel=photons_per_pixel,
        rng=np.random.default_rng(seed),
        read_noise_electrons=read_noise_electrons,
        bin_size=15,
        return_intermediates=True,
    )

    np.testing.assert_allclose(result["binned_image"], expected_binned)
    np.testing.assert_allclose(result["noisy_counts"], expected_counts)
    np.testing.assert_allclose(result["noisy_image"], expected_image)


def test_simulate_noisy_camera_image_can_return_counts_for_binned_input() -> None:
    binned = np.array([[0.5, 1.0], [1.5, 2.0]])
    photons_per_pixel = 80.0
    read_noise_electrons = 1.25
    seed = 19

    expected_counts = add_camera_noise(
        binned,
        photons_per_pixel,
        np.random.default_rng(seed),
        read_noise_electrons,
    )
    result = simulate_noisy_camera_image(
        binned,
        photons_per_pixel=photons_per_pixel,
        rng=np.random.default_rng(seed),
        read_noise_electrons=read_noise_electrons,
        input_is_binned=True,
        normalize=False,
    )

    np.testing.assert_allclose(result, expected_counts)
