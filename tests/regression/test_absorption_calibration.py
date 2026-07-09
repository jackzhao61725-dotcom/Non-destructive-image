from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import (
    compute_optical_density,
    estimate_cloud_moments,
    extract_absorption_observables,
    integrate_optical_density,
)


def test_compute_optical_density_matches_known_atom_probe_dark_images() -> None:
    atom = np.full((2, 2), 50.0)
    probe = np.full((2, 2), 100.0)
    dark = np.full((2, 2), 10.0)

    od = compute_optical_density(atom, probe, dark)

    expected = -np.log((50.0 - 10.0) / (100.0 - 10.0))
    np.testing.assert_allclose(od, expected)


def test_compute_optical_density_clipping_avoids_log_of_zero() -> None:
    atom = np.array([[0.0, 1.0]])
    probe = np.array([[0.0, 2.0]])

    od = compute_optical_density(atom, probe, epsilon=1e-6)

    assert od.shape == atom.shape
    assert np.all(np.isfinite(od))


def test_integrated_optical_density_is_finite_and_area_scaled() -> None:
    od = np.array([[1.0, 2.0], [3.0, 4.0]])

    assert integrate_optical_density(od) == pytest.approx(10.0)
    assert integrate_optical_density(od, pixel_area=0.25) == pytest.approx(2.5)


def test_cloud_moments_match_simple_symmetric_od_map() -> None:
    od = np.array(
        [
            [1.0, 2.0, 1.0],
            [2.0, 4.0, 2.0],
            [1.0, 2.0, 1.0],
        ]
    )
    x = np.array([-1.0, 0.0, 1.0])
    y = np.array([-2.0, 0.0, 2.0])

    moments = estimate_cloud_moments(od, x=x, y=y)

    assert moments["peak_od"] == pytest.approx(4.0)
    assert moments["integrated_od"] == pytest.approx(16.0)
    assert moments["centre_x"] == pytest.approx(0.0)
    assert moments["centre_y"] == pytest.approx(0.0)
    assert moments["width_x"] == pytest.approx(np.sqrt(0.5))
    assert moments["width_y"] == pytest.approx(np.sqrt(2.0))


def test_extract_absorption_observables_is_deterministic() -> None:
    probe = np.full((3, 3), 100.0)
    od_target = np.array(
        [
            [0.1, 0.2, 0.1],
            [0.2, 0.4, 0.2],
            [0.1, 0.2, 0.1],
        ]
    )
    atom = probe * np.exp(-od_target)

    first = extract_absorption_observables(atom, probe)
    second = extract_absorption_observables(atom, probe)

    assert first == second
    assert set(first) == {"peak_od", "integrated_od", "centre_x", "centre_y", "width_x", "width_y"}


def test_absorption_helpers_reject_invalid_shapes() -> None:
    atom = np.ones((2, 2))
    probe = np.ones((2, 3))

    with pytest.raises(ValueError, match="probe_image must have shape"):
        compute_optical_density(atom, probe)


def test_extract_absorption_observables_applies_pixel_area_to_integrated_od() -> None:
    probe = np.full((2, 2), 100.0)
    atom = probe * np.exp(-0.5)

    observables = extract_absorption_observables(atom, probe, pixel_area=0.2)

    assert observables["integrated_od"] == pytest.approx(4 * 0.5 * 0.2)
