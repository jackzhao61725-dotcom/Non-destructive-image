from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from non_destructive_image import simulate_faraday_image
from non_destructive_image.camera import bin_to_camera_pixels
from non_destructive_image.reconstruction import (
    DarkFieldFaradayMeasurement,
    DetectorContract,
    DualPortFaradayMeasurement,
    FaradayResponseContract,
    FitOptions,
    ReconstructionGrid,
    SmoothTFBounds,
    SmoothTFParameters,
    assert_response_scale_degeneracy,
    estimate_dual_port_initial_parameters,
    fit_smooth_tf,
)
from non_destructive_image.reconstruction.parameters import from_internal, to_internal


@pytest.fixture(scope="module")
def reconstruction_context() -> dict[str, object]:
    size = 64
    field_of_view_m = 80e-6
    spacing_m = field_of_view_m / size
    axis_m = (np.arange(size) - size / 2) * spacing_m
    y_grid_m, z_grid_m = np.meshgrid(axis_m, axis_m)
    frequency = np.fft.fftfreq(size, d=spacing_m)
    fy, fz = np.meshgrid(frequency, frequency)
    pupil = ((fy**2 + fz**2) <= (2.2e5) ** 2).astype(float)
    bin_size = 4
    camera_y_um = bin_to_camera_pixels(y_grid_m * 1e6, bin_size)[0, :]
    camera_z_um = bin_to_camera_pixels(z_grid_m * 1e6, bin_size)[:, 0]
    camera_y, camera_z = np.meshgrid(camera_y_um, camera_z_um)
    roi_mask = (np.abs(camera_y) <= 27.0) & (np.abs(camera_z) <= 16.0)
    grid = ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=bin_size,
        roi_mask=roi_mask,
    )
    detector = DetectorContract(
        photoelectrons_per_i0_pixel=900.0,
        read_noise_electrons_per_pixel_per_readout=3.0,
    )
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=4.0e-16,
        kappa_f=1.0,
    )
    truth = SmoothTFParameters(
        column_density_peak_m2=5.0e14,
        y0_um=1.0,
        z0_um=-0.5,
        radius_y_um=18.0,
        radius_z_um=5.0,
    )
    bounds = SmoothTFBounds.from_mapping(
        {
            "column_density_peak_m2": [1.0e13, 5.0e15],
            "y0_um": [-10.0, 10.0],
            "z0_um": [-8.0, 8.0],
            "radius_y_um": [5.0, 50.0],
            "radius_z_um": [0.5, 15.0],
        }
    )
    return {
        "grid": grid,
        "detector": detector,
        "response": response,
        "truth": truth,
        "bounds": bounds,
        "dual": DualPortFaradayMeasurement(grid=grid, detector=detector, response=response),
        "dark": DarkFieldFaradayMeasurement(grid=grid, detector=detector, response=response),
    }


def test_reconstruction_contract_uses_absolute_not_truth_relative_bounds() -> None:
    config_path = Path(__file__).parents[2] / "configs" / "reconstruction_v3.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    bounds = SmoothTFBounds.from_mapping(config["parameter_bounds"])
    assert bounds.lower.column_density_peak_m2 == pytest.approx(1.0e13)
    assert bounds.upper.radius_y_um == pytest.approx(60.0)
    assert not any("truth" in key.lower() for key in config["parameter_bounds"])


def test_truth_relative_bound_keys_are_rejected() -> None:
    with pytest.raises(ValueError, match="must not depend on simulated truth"):
        SmoothTFBounds.from_mapping(
            {
                "column_density_peak_fraction_of_truth": [0.2, 2.5],
                "column_density_peak_m2": [1.0e13, 5.0e15],
                "y0_um": [-10.0, 10.0],
                "z0_um": [-8.0, 8.0],
                "radius_y_um": [5.0, 50.0],
                "radius_z_um": [0.5, 15.0],
            }
        )


def test_raw_faraday_measurements_match_existing_forward_model(
    reconstruction_context: dict[str, object],
) -> None:
    dual = reconstruction_context["dual"]
    dark = reconstruction_context["dark"]
    truth = reconstruction_context["truth"]
    assert isinstance(dual, DualPortFaradayMeasurement)
    assert isinstance(dark, DarkFieldFaradayMeasurement)
    assert isinstance(truth, SmoothTFParameters)
    theta = dual.response.rotation_per_column_density_rad_m2 * dual.column_density(truth)
    complete = simulate_faraday_image(theta, dual.grid.pupil)
    expected_h = bin_to_camera_pixels(complete["dual_port_v_intensity"], dual.grid.bin_size)
    expected_v = bin_to_camera_pixels(complete["dual_port_u_intensity"], dual.grid.bin_size)
    actual_h, actual_v = dual.normalised_channels(truth)
    np.testing.assert_allclose(actual_h, expected_h, rtol=0.0, atol=5e-13)
    np.testing.assert_allclose(actual_v, expected_v, rtol=0.0, atol=5e-13)
    (actual_dark,) = dark.normalised_channels(truth)
    expected_dark = bin_to_camera_pixels(complete["dark_field_intensity"], dark.grid.bin_size)
    np.testing.assert_allclose(actual_dark, expected_dark, rtol=0.0, atol=5e-13)


def test_dual_port_initialisation_tracks_a_negative_faraday_response(
    reconstruction_context: dict[str, object],
) -> None:
    grid = reconstruction_context["grid"]
    detector = reconstruction_context["detector"]
    truth = reconstruction_context["truth"]
    bounds = reconstruction_context["bounds"]
    assert isinstance(grid, ReconstructionGrid)
    assert isinstance(detector, DetectorContract)
    assert isinstance(truth, SmoothTFParameters)
    assert isinstance(bounds, SmoothTFBounds)
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=4.0e-16,
        kappa_f=-45.0 / 91.0,
    )
    measurement = DualPortFaradayMeasurement(
        grid=grid,
        detector=detector,
        response=response,
    )
    observed_h, observed_v = measurement.expected_channels(truth)

    initial = estimate_dual_port_initial_parameters(
        measurement,
        observed_h,
        observed_v,
        bounds,
    )
    positive_measurement = DualPortFaradayMeasurement(
        grid=grid,
        detector=detector,
        response=FaradayResponseContract(
            phase_per_column_density_rad_m2=4.0e-16,
            kappa_f=45.0 / 91.0,
        ),
    )
    positive_h, positive_v = positive_measurement.expected_channels(truth)
    positive_initial = estimate_dual_port_initial_parameters(
        positive_measurement,
        positive_h,
        positive_v,
        bounds,
    )

    assert initial.y0_um == pytest.approx(truth.y0_um, abs=1.5)
    assert initial.z0_um == pytest.approx(truth.z0_um, abs=1.5)
    assert initial.column_density_peak_m2 > bounds.lower.column_density_peak_m2
    np.testing.assert_allclose(
        to_internal(initial),
        to_internal(positive_initial),
        rtol=0.0,
        atol=1e-12,
    )


@pytest.mark.parametrize("measurement_key", ["dual", "dark"])
def test_analytic_jacobian_matches_central_difference(
    reconstruction_context: dict[str, object],
    measurement_key: str,
) -> None:
    measurement = reconstruction_context[measurement_key]
    truth = reconstruction_context["truth"]
    assert isinstance(
        measurement,
        (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
    )
    assert isinstance(truth, SmoothTFParameters)
    internal = to_internal(truth)
    _, analytic = measurement.expected_vector_and_jacobian_internal(internal)
    steps = np.asarray([2e-6, 2e-5, 2e-5, 2e-6, 2e-6])
    for index, step in enumerate(steps):
        plus = internal.copy()
        minus = internal.copy()
        plus[index] += step
        minus[index] -= step
        plus_vector = measurement.flatten_observed(
            *measurement.expected_channels(from_internal(plus))
        )
        minus_vector = measurement.flatten_observed(
            *measurement.expected_channels(from_internal(minus))
        )
        numerical = (plus_vector - minus_vector) / (2.0 * step)
        relative_error = np.linalg.norm(numerical - analytic[:, index]) / max(
            np.linalg.norm(numerical),
            np.finfo(float).eps,
        )
        assert relative_error < 3e-6


@pytest.mark.parametrize("measurement_key", ["dual", "dark"])
def test_density_and_kappa_scale_are_structurally_degenerate(
    reconstruction_context: dict[str, object],
    measurement_key: str,
) -> None:
    measurement = reconstruction_context[measurement_key]
    truth = reconstruction_context["truth"]
    assert isinstance(
        measurement,
        (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
    )
    assert isinstance(truth, SmoothTFParameters)
    assert_response_scale_degeneracy(measurement, truth, scale=0.5)


@pytest.mark.parametrize("measurement_key", ["dual", "dark"])
def test_noiseless_truth_near_fit_recovers_the_smooth_object(
    reconstruction_context: dict[str, object],
    measurement_key: str,
) -> None:
    measurement = reconstruction_context[measurement_key]
    truth = reconstruction_context["truth"]
    bounds = reconstruction_context["bounds"]
    assert isinstance(
        measurement,
        (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
    )
    assert isinstance(truth, SmoothTFParameters)
    assert isinstance(bounds, SmoothTFBounds)
    start = SmoothTFParameters(
        column_density_peak_m2=truth.column_density_peak_m2 * 1.03,
        y0_um=truth.y0_um + 0.1,
        z0_um=truth.z0_um - 0.1,
        radius_y_um=truth.radius_y_um * 0.98,
        radius_z_um=truth.radius_z_um * 1.03,
    )
    result = fit_smooth_tf(
        measurement,
        measurement.expected_channels(truth),
        starts=[start],
        bounds=bounds,
        options=FitOptions(irls_iterations=2, max_nfev=80),
    )
    assert result.diagnostics.success
    assert result.diagnostics.jacobian_rank == 5
    np.testing.assert_allclose(
        to_internal(result.parameters),
        to_internal(truth),
        rtol=0.0,
        atol=2e-5,
    )
