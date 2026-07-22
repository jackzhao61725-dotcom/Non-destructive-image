"""End-to-end tests for provenance-bearing experimental reconstruction."""

from __future__ import annotations

import hashlib
from dataclasses import fields

import numpy as np
import pytest

from non_destructive_image.reconstruction.contracts import (
    DetectorContract,
    FaradayResponseContract,
    ReconstructionGrid,
)
from non_destructive_image.reconstruction.density_fit import DensityFitOptions
from non_destructive_image.reconstruction.evidence import ResponseScaleDeclaration
from non_destructive_image.reconstruction.experimental_adapter import CameraCalibration
from non_destructive_image.reconstruction.experimental_reconstruction import (
    BootstrapRequest,
    reconstruct_experimental_exposure,
)
from non_destructive_image.reconstruction.measurements import (
    DualPortFaradayMeasurement,
)
from non_destructive_image.reconstruction.object_models import (
    NonnegativeBilinearDensityModel,
)
from non_destructive_image.reconstruction.observables import (
    ObservableIntegrationSupport,
    extract_density_observables,
)
from non_destructive_image.reconstruction.observations import (
    ExperimentContext,
    NonDestructiveExposure,
    NonDestructiveSequence,
    ProbePulse,
    RawCameraFrame,
    UncertainQuantity,
)


def _sha256(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


@pytest.fixture(scope="module")
def reconstruction_case() -> dict[str, object]:
    coordinate = np.linspace(-4e-6, 4e-6, 8)
    y_grid, z_grid = np.meshgrid(coordinate, coordinate)
    grid = ReconstructionGrid.from_arrays(
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        pupil=np.ones((8, 8), dtype=complex),
        bin_size=2,
        roi_mask=np.ones((4, 4), dtype=bool),
    )
    measurement = DualPortFaradayMeasurement(
        grid=grid,
        detector=DetectorContract(1000.0, 3.0),
        response=FaradayResponseContract(1.0e-14, 1.0),
    )
    model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        knot_y_um=(-4.0, 0.0, 4.0),
        knot_z_um=(-4.0, 0.0, 4.0),
        coefficient_scale_m2=1.0e13,
    )
    observable_support = ObservableIntegrationSupport(
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        support_mask=(np.abs(y_grid) <= 2.0e-6) & (np.abs(z_grid) <= 2.0e-6),
    )
    truth_coefficients = np.asarray(
        [
            0.22,
            0.30,
            0.20,
            0.34,
            0.58,
            0.31,
            0.18,
            0.27,
            0.19,
        ]
    )
    expected_h, expected_v = measurement.expected_channels_from_density(
        model.column_density(truth_coefficients)
    )
    quantity = UncertainQuantity(1.0, 0.1, "1", source="test declaration")
    context = ExperimentContext(
        context_id="context-electron-test",
        condensate_quantities={"population_scale": quantity},
        preparation_quantities={"preparation_scale": quantity},
        instrument_quantities={"instrument_scale": quantity},
    )
    pulse = ProbePulse(
        pulse_id="pulse-0",
        wavelength_m=401e-9,
        detuning_hz=-1.5e9,
        power_w=1.0e-3,
        duration_s=40e-6,
        polarisation_label="linear",
    )
    exposure = NonDestructiveExposure(
        exposure_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        readout="dual_port",
        # Deliberately reverse insertion order: the adapter owns channel order.
        raw_channels={
            "V": RawCameraFrame(expected_v, "electron", "camera-1"),
            "H": RawCameraFrame(expected_h, "electron", "camera-1"),
        },
    )
    sequence = NonDestructiveSequence(
        sequence_id="sequence-electron-test",
        context=context,
        exposures=(exposure,),
    )
    return {
        "measurement": measurement,
        "model": model,
        "observable_support": observable_support,
        "truth_coefficients": truth_coefficients,
        "expected_h": expected_h,
        "expected_v": expected_v,
        "context": context,
        "pulse": pulse,
        "sequence": sequence,
    }


def test_electron_sequence_returns_uncalibrated_evidence_bundle_without_truth(
    reconstruction_case: dict[str, object],
) -> None:
    measurement = reconstruction_case["measurement"]
    model = reconstruction_case["model"]
    observable_support = reconstruction_case["observable_support"]
    truth_coefficients = reconstruction_case["truth_coefficients"]
    sequence = reconstruction_case["sequence"]
    expected_h = reconstruction_case["expected_h"]
    expected_v = reconstruction_case["expected_v"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(observable_support, ObservableIntegrationSupport)
    assert isinstance(sequence, NonDestructiveSequence)

    result = reconstruct_experimental_exposure(
        sequence,
        0,
        measurement=measurement,
        camera_calibrations={
            "camera-1": CameraCalibration(
                camera_id="camera-1",
                electrons_per_adu=2.0,
                read_noise_electrons=3.0,
                bias_adu=10.0,
            )
        },
        model=model,
        observable_integration_support=observable_support,
        initial_coefficients=truth_coefficients,
        coefficient_upper=1.5,
        regularisation=None,
        response_scale=ResponseScaleDeclaration(
            kappa_f=1.0,
            status="illustrative_uncalibrated",
            source="declared test operator",
        ),
        pipeline_fingerprint="basis3x3-fit-v1",
        condition_fingerprint="dual-electron-F90-v1",
        fit_options=DensityFitOptions(irls_iterations=1, max_nfev=40),
        bootstrap_request=BootstrapRequest(draws=2, seed=20260722),
    )

    assert result.context_id == "context-electron-test"
    assert result.sequence_id == "sequence-electron-test"
    assert result.channel_names == ("H", "V")
    assert tuple(name for name, _ in result.input_channel_sha256) == ("H", "V")
    assert dict(result.input_channel_sha256) == {
        "H": _sha256(expected_h),
        "V": _sha256(expected_v),
    }
    assert result.response_scale.status == "illustrative_uncalibrated"
    assert result.absolute_density_status == "conditional_on_assumed_kappa_f"
    assert result.fit_diagnostics.success
    assert result.latent_fit is None
    assert result.local_credibility.density_standard_uncertainty_m2 is None
    assert result.observables.integration_support.is_identical_to(observable_support)
    assert result.null_evidence.target_origin == "experimental_observation"
    assert result.null_evidence.evidence_level == "model_only"
    assert result.bootstrap.status == "complete"
    assert result.bootstrap.result is not None
    assert result.bootstrap.message is not None
    assert result.bootstrap.result.coefficient_samples is None
    assert result.bootstrap.result.density_mean_m2 is None
    assert result.bootstrap.result.density_standard_uncertainty_m2 is None
    assert result.bootstrap.result.density_interval_lower_m2 is None
    assert result.bootstrap.result.density_interval_upper_m2 is None
    assert result.bootstrap.result.feature_intervals == {}
    assert result.bootstrap.result.observables is not None
    assert result.bootstrap.result.observables.samples.shape == (2, 4)
    bootstrap_point_support = (
        result.bootstrap.result.observables.point_estimate.integration_support
    )
    assert (
        bootstrap_point_support.is_identical_to(observable_support)
    )
    assert result.stability is None
    assert any("No truth density" in item for item in result.interpretation_limits)
    assert any(
        "conditional on the assumed kappa_F" in item
        for item in result.interpretation_limits
    )

    result_fields = {definition.name for definition in fields(type(result))}
    assert not any("truth" in name.lower() for name in result_fields)
    assert not hasattr(result, "truth_density")
    assert not hasattr(result, "reference_density")
    assert not hasattr(result, "fit")
    assert not hasattr(result, "features")

    diagnostic_result = reconstruct_experimental_exposure(
        sequence,
        0,
        measurement=measurement,
        camera_calibrations={
            "camera-1": CameraCalibration(
                camera_id="camera-1",
                electrons_per_adu=2.0,
                read_noise_electrons=3.0,
                bias_adu=10.0,
            )
        },
        model=model,
        observable_integration_support=observable_support,
        initial_coefficients=truth_coefficients,
        coefficient_upper=1.5,
        regularisation=None,
        response_scale=ResponseScaleDeclaration(
            kappa_f=1.0,
            status="illustrative_uncalibrated",
            source="declared test operator",
        ),
        pipeline_fingerprint="basis3x3-fit-v1",
        condition_fingerprint="dual-electron-F90-v1",
        fit_options=DensityFitOptions(irls_iterations=1, max_nfev=40),
        retain_latent_artifacts=True,
    )
    assert diagnostic_result.latent_fit is not None
    assert (
        diagnostic_result.local_credibility.density_standard_uncertainty_m2
        is not None
    )
    expected_observables = extract_density_observables(
        diagnostic_result.latent_fit.column_density_m2,
        observable_support,
    )
    np.testing.assert_allclose(
        [
            result.observables.integrated_response,
            result.observables.centroid_y_m,
            result.observables.centroid_z_m,
            result.observables.major_rms_width_m,
        ],
        [
            expected_observables.integrated_response,
            expected_observables.centroid_y_m,
            expected_observables.centroid_z_m,
            expected_observables.major_rms_width_m,
        ],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        [
            diagnostic_result.observables.integrated_response,
            diagnostic_result.observables.centroid_y_m,
            diagnostic_result.observables.centroid_z_m,
            diagnostic_result.observables.major_rms_width_m,
        ],
        [
            result.observables.integrated_response,
            result.observables.centroid_y_m,
            result.observables.centroid_z_m,
            result.observables.major_rms_width_m,
        ],
        rtol=0.0,
        atol=0.0,
    )
    full_support = ObservableIntegrationSupport(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
    )
    full_observables = extract_density_observables(
        diagnostic_result.latent_fit.column_density_m2,
        full_support,
    )
    assert result.observables.integrated_response < full_observables.integrated_response


def test_observable_support_coordinates_and_latent_extent_are_strict(
    reconstruction_case: dict[str, object],
) -> None:
    measurement = reconstruction_case["measurement"]
    model = reconstruction_case["model"]
    observable_support = reconstruction_case["observable_support"]
    truth_coefficients = reconstruction_case["truth_coefficients"]
    sequence = reconstruction_case["sequence"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(observable_support, ObservableIntegrationSupport)
    assert isinstance(sequence, NonDestructiveSequence)

    common = dict(
        measurement=measurement,
        camera_calibrations={},
        model=model,
        initial_coefficients=truth_coefficients,
        coefficient_upper=1.5,
        regularisation=None,
        response_scale=ResponseScaleDeclaration(
            kappa_f=1.0,
            status="illustrative_uncalibrated",
            source="declared test operator",
        ),
        pipeline_fingerprint="basis3x3-fit-v1",
        condition_fingerprint="dual-electron-F90-v1",
        fit_options=DensityFitOptions(irls_iterations=1, max_nfev=40),
    )
    shifted_support = ObservableIntegrationSupport(
        y_grid_m=observable_support.y_grid_m + 1.0e-12,
        z_grid_m=observable_support.z_grid_m,
        support_mask=observable_support.support_mask,
        cell_area_m2=observable_support.cell_area_m2,
    )
    with pytest.raises(ValueError, match="coordinates must exactly match"):
        reconstruct_experimental_exposure(
            sequence,
            0,
            observable_integration_support=shifted_support,
            **common,
        )

    restricted_model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
        knot_y_um=model.knot_y_um,
        knot_z_um=model.knot_z_um,
        coefficient_scale_m2=model.coefficient_scale_m2,
        support_mask=observable_support.support_mask,
    )
    full_support = ObservableIntegrationSupport(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
    )
    with pytest.raises(ValueError, match="cannot extend beyond latent model support"):
        reconstruct_experimental_exposure(
            sequence,
            0,
            observable_integration_support=full_support,
            **{**common, "model": restricted_model},
        )
    with pytest.raises(ValueError, match="stability density maps require"):
        reconstruct_experimental_exposure(
            sequence,
            0,
            observable_integration_support=observable_support,
            stability_variants={},
            **common,
        )


def test_end_to_end_adu_input_requires_exact_camera_calibration(
    reconstruction_case: dict[str, object],
) -> None:
    measurement = reconstruction_case["measurement"]
    model = reconstruction_case["model"]
    observable_support = reconstruction_case["observable_support"]
    truth_coefficients = reconstruction_case["truth_coefficients"]
    context = reconstruction_case["context"]
    pulse = reconstruction_case["pulse"]
    expected_h = reconstruction_case["expected_h"]
    expected_v = reconstruction_case["expected_v"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(observable_support, ObservableIntegrationSupport)
    assert isinstance(context, ExperimentContext)
    assert isinstance(pulse, ProbePulse)

    adu_sequence = NonDestructiveSequence(
        sequence_id="sequence-adu-test",
        context=context,
        exposures=(
            NonDestructiveExposure(
                exposure_index=0,
                timestamp_s=0.0,
                pulse=pulse,
                readout="dual_port",
                raw_channels={
                    "H": RawCameraFrame(expected_h / 2.0 + 10.0, "adu", "camera-1"),
                    "V": RawCameraFrame(expected_v / 2.0 + 10.0, "adu", "camera-1"),
                },
            ),
        ),
    )
    common = dict(
        measurement=measurement,
        model=model,
        observable_integration_support=observable_support,
        initial_coefficients=truth_coefficients,
        coefficient_upper=1.5,
        regularisation=None,
        response_scale=ResponseScaleDeclaration(
            kappa_f=1.0,
            status="illustrative_uncalibrated",
            source="declared test operator",
        ),
        pipeline_fingerprint="basis3x3-fit-v1",
        condition_fingerprint="dual-adu-F90-v1",
        fit_options=DensityFitOptions(irls_iterations=1, max_nfev=40),
    )

    with pytest.raises(ValueError, match="missing=.*camera-1"):
        reconstruct_experimental_exposure(
            adu_sequence,
            0,
            camera_calibrations={},
            **common,
        )
    with pytest.raises(ValueError, match="read noise .* incompatible"):
        reconstruct_experimental_exposure(
            adu_sequence,
            0,
            camera_calibrations={
                "camera-1": CameraCalibration("camera-1", 2.0, 2.5, 10.0)
            },
            **common,
        )
