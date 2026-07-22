from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.reconstruction import (
    DarkFieldFaradayMeasurement,
    DetectorContract,
    DualPortFaradayMeasurement,
    FaradayResponseContract,
    NonnegativeBilinearDensityModel,
    build_physical_camera_aligned_reduced_grid,
    build_uniform_physical_camera_grid,
)


def test_orca_m10_camera_grid_preserves_constants_and_physical_coordinates() -> None:
    grid = build_uniform_physical_camera_grid(
        ngrid=1024,
        field_of_view_m=100e-6,
        camera_pixel_size_m=0.650e-6,
        camera_output_shape=(153, 153),
        numerical_aperture=0.080,
        wavelength_m=401e-9,
    )

    sampled = grid.camera_average(np.ones((1024, 1024), dtype=float))

    assert grid.sampling_mode == "physical_pixel"
    assert grid.camera_shape == (153, 153)
    np.testing.assert_allclose(sampled, 1.0, rtol=0.0, atol=1e-13)
    assert grid.camera_y_um[76] == 0.0
    assert grid.camera_z_um[76] == 0.0
    np.testing.assert_allclose(np.diff(grid.camera_y_um), 0.650, atol=1e-14)
    np.testing.assert_allclose(np.diff(grid.camera_z_um), 0.650, atol=1e-14)
    np.testing.assert_allclose(grid.camera_y_um[[0, -1]], [-49.4, 49.4])


def test_physical_camera_reduced_grid_centres_match_orca_sampling() -> None:
    reduced, provenance = build_physical_camera_aligned_reduced_grid(
        canonical_ngrid=1024,
        canonical_field_of_view_m=100e-6,
        camera_pixel_size_m=0.650e-6,
        camera_output_shape=(153, 153),
        reduced_ngrid=306,
        numerical_aperture=0.080,
        wavelength_m=401e-9,
    )

    assert reduced.camera_shape == (153, 153)
    assert provenance.reduced_ngrid == 306
    assert provenance.reduced_field_of_view_m == 100e-6
    assert provenance.maximum_camera_coordinate_mismatch_m < 1e-15
    np.testing.assert_allclose(np.diff(reduced.camera_y_um), 0.650, atol=1e-14)


@pytest.mark.parametrize(
    "measurement_type",
    [DualPortFaradayMeasurement, DarkFieldFaradayMeasurement],
)
def test_physical_camera_stack_and_jacobian_use_the_same_sampler(
    measurement_type: type[DualPortFaradayMeasurement]
    | type[DarkFieldFaradayMeasurement],
) -> None:
    grid = build_uniform_physical_camera_grid(
        ngrid=64,
        field_of_view_m=80e-6,
        camera_pixel_size_m=4.7e-6,
        camera_output_shape=(15, 15),
        numerical_aperture=0.12,
        wavelength_m=401e-9,
    )
    rng = np.random.default_rng(20260721)
    object_stack = rng.normal(size=(3, 64, 64))
    sampled_stack = grid.camera_average_stack(object_stack)
    expected_stack = np.stack(
        [grid.camera_average(image) for image in object_stack],
        axis=0,
    )
    np.testing.assert_allclose(sampled_stack, expected_stack, rtol=0.0, atol=1e-14)

    measurement = measurement_type(
        grid=grid,
        detector=DetectorContract(
            photoelectrons_per_i0_pixel=220.58087277528466,
            read_noise_electrons_per_pixel_per_readout=1.4,
        ),
        response=FaradayResponseContract(
            phase_per_column_density_rad_m2=3.77498270624925e-16,
            kappa_f=1.0,
        ),
    )
    model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=[-18.0, -9.0, 0.0, 9.0, 18.0],
        knot_z_um=[-6.0, 0.0, 6.0],
        coefficient_scale_m2=5.0e14,
    )
    coefficients = np.asarray(
        [
            [0.02, 0.08, 0.12, 0.05, 0.01],
            [0.05, 0.90, 0.18, 0.72, 0.04],
            [0.01, 0.06, 0.14, 0.09, 0.02],
        ],
        dtype=float,
    ).ravel()
    _, analytic = measurement.expected_vector_and_jacobian_model(model, coefficients)
    index = 7
    step = 2e-6
    plus = coefficients.copy()
    minus = coefficients.copy()
    plus[index] += step
    minus[index] -= step
    numerical = (
        measurement.expected_vector_from_density(model.column_density(plus))
        - measurement.expected_vector_from_density(model.column_density(minus))
    ) / (2.0 * step)
    relative_error = np.linalg.norm(numerical - analytic[:, index]) / max(
        np.linalg.norm(numerical),
        np.finfo(float).eps,
    )
    assert relative_error < 5e-6
