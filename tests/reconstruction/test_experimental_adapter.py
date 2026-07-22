"""Tests for strict camera-to-reconstruction experimental adapters."""

from __future__ import annotations

from dataclasses import fields

import numpy as np
import pytest

from non_destructive_image.reconstruction.contracts import (
    DetectorContract,
    FaradayResponseContract,
    ReconstructionGrid,
)
from non_destructive_image.reconstruction.experimental_adapter import (
    CameraCalibration,
    bind_non_destructive_exposure,
    build_rai_transmission,
    calibrate_camera_frame,
)
from non_destructive_image.reconstruction.measurements import (
    DarkFieldFaradayMeasurement,
    DualPortFaradayMeasurement,
)
from non_destructive_image.reconstruction.observations import (
    ExperimentContext,
    NonDestructiveExposure,
    PreNormalisedCameraImage,
    ProbePulse,
    RawCameraFrame,
    ResonantAbsorptionObservation,
    UncertainQuantity,
)


@pytest.fixture
def context() -> ExperimentContext:
    quantity = UncertainQuantity(1.0, 0.1, "1")
    return ExperimentContext(
        context_id="adapter-test",
        condensate_quantities={"population_scale": quantity},
        preparation_quantities={"preparation_scale": quantity},
        instrument_quantities={"instrument_scale": quantity},
    )


@pytest.fixture
def pulse() -> ProbePulse:
    return ProbePulse(
        pulse_id="pulse",
        wavelength_m=401e-9,
        detuning_hz=-1.5e9,
        power_w=1.0e-3,
        duration_s=40e-6,
        polarisation_label="linear",
    )


@pytest.fixture
def grid() -> ReconstructionGrid:
    coordinate = np.linspace(-4e-6, 4e-6, 8)
    y_grid, z_grid = np.meshgrid(coordinate, coordinate)
    return ReconstructionGrid.from_arrays(
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        pupil=np.ones((8, 8), dtype=complex),
        bin_size=2,
        roi_mask=np.asarray(
            [
                [False, True, True, False],
                [True, True, True, True],
                [True, True, True, True],
                [False, True, True, False],
            ]
        ),
    )


@pytest.fixture
def dual_measurement(grid: ReconstructionGrid) -> DualPortFaradayMeasurement:
    return DualPortFaradayMeasurement(
        grid=grid,
        detector=DetectorContract(1000.0, 3.0),
        response=FaradayResponseContract(1e-14, 1.0),
    )


@pytest.fixture
def dark_measurement(grid: ReconstructionGrid) -> DarkFieldFaradayMeasurement:
    return DarkFieldFaradayMeasurement(
        grid=grid,
        detector=DetectorContract(1000.0, 3.0),
        response=FaradayResponseContract(1e-14, 1.0),
    )


def frame(
    values: float | np.ndarray,
    *,
    unit: str = "electron",
    camera_id: str = "camera-1",
    shape: tuple[int, int] = (4, 4),
) -> RawCameraFrame:
    array = np.full(shape, values, dtype=float) if np.isscalar(values) else np.asarray(values)
    return RawCameraFrame(array, unit, camera_id)  # type: ignore[arg-type]


def calibration(
    *,
    camera_id: str = "camera-1",
    gain: float = 2.0,
    read_noise: float = 3.0,
    bias: float | np.ndarray = 10.0,
) -> CameraCalibration:
    return CameraCalibration(camera_id, gain, read_noise, bias)


def rai_raw(
    context: ExperimentContext,
    pulse: ProbePulse,
    *,
    atom: RawCameraFrame,
    reference: RawCameraFrame,
    dark: RawCameraFrame,
) -> ResonantAbsorptionObservation:
    return ResonantAbsorptionObservation(
        context=context,
        observation_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        atom_frame=atom,
        reference_frame=reference,
        dark_frame=dark,
    )


def test_camera_calibration_converts_adu_with_scalar_bias() -> None:
    converted = calibrate_camera_frame(frame(14.0, unit="adu"), calibration())
    np.testing.assert_allclose(converted.values_electrons, 8.0)
    assert converted.source_unit == "adu"
    assert not converted.values_electrons.flags.writeable
    assert "electrons-per-ADU" in converted.assumptions[0]


def test_camera_calibration_supports_bias_map_and_requires_matching_shape() -> None:
    values = np.arange(16.0).reshape(4, 4) + 20.0
    bias = np.arange(16.0).reshape(4, 4)
    converted = calibrate_camera_frame(
        frame(values, unit="adu"),
        calibration(gain=0.5, bias=bias),
    )
    np.testing.assert_allclose(converted.values_electrons, 10.0)
    with pytest.raises(ValueError, match="bias map has shape"):
        calibrate_camera_frame(
            frame(values, unit="adu"),
            calibration(bias=np.ones((3, 4))),
        )


def test_electron_input_is_not_recalibrated() -> None:
    values = np.arange(16.0).reshape(4, 4) - 3.0
    converted = calibrate_camera_frame(
        frame(values, unit="electron"),
        calibration(gain=100.0, bias=np.ones((2, 2)) * 999.0),
    )
    np.testing.assert_array_equal(converted.values_electrons, values)
    assert "no ADU conversion" in converted.assumptions[0]
    assert "cannot independently verify" in converted.warnings[0]


def test_adu_requires_matching_valid_calibration() -> None:
    with pytest.raises(ValueError, match="requires an explicit"):
        calibrate_camera_frame(frame(10.0, unit="adu"), None)
    with pytest.raises(ValueError, match="does not match calibration"):
        calibrate_camera_frame(
            frame(10.0, unit="adu"),
            calibration(camera_id="camera-2"),
        )
    with pytest.raises(ValueError, match="positive"):
        calibration(gain=0.0)
    with pytest.raises(ValueError, match="cannot be negative"):
        calibration(read_noise=-1.0)


def test_raw_rai_triplet_builds_transmission_not_density(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    observation = rai_raw(
        context,
        pulse,
        atom=frame(50.0, unit="adu"),
        reference=frame(90.0, unit="adu"),
        dark=frame(10.0, unit="adu"),
    )
    bundle = build_rai_transmission(observation, calibration())
    np.testing.assert_allclose(bundle.transmission, 0.5)
    np.testing.assert_allclose(bundle.atom_electrons, 80.0)
    np.testing.assert_allclose(bundle.reference_electrons, 160.0)
    np.testing.assert_allclose(bundle.dark_electrons, 0.0)
    assert bundle.representation == "raw_triplet"
    assert np.all(bundle.valid_mask)
    assert any("No Beer--Lambert" in statement for statement in bundle.assumptions)
    field_names = {definition.name for definition in fields(type(bundle))}
    assert "density" not in field_names
    assert "optical_depth" not in field_names


def test_electron_rai_triplet_needs_no_adu_calibration(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    observation = rai_raw(
        context,
        pulse,
        atom=frame(55.0),
        reference=frame(105.0),
        dark=frame(5.0),
    )
    bundle = build_rai_transmission(observation)
    np.testing.assert_allclose(bundle.transmission, 0.5)
    assert any("upstream electron calibration" in item for item in bundle.warnings)


def test_rai_rejects_nonpositive_reference_minus_dark(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    observation = rai_raw(
        context,
        pulse,
        atom=frame(4.0),
        reference=frame(5.0),
        dark=frame(5.0),
    )
    with pytest.raises(ValueError, match="reference-minus-dark must be positive"):
        build_rai_transmission(observation)


def test_pre_normalised_rai_requires_explicit_transmission_semantics(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    image = PreNormalisedCameraImage(
        np.asarray([[0.8, 1.1], [-0.1, 0.9]]),
        quantity_name="normalised transmission",
        normalisation_description="registered (atom-dark)/(reference-dark)",
    )
    observation = ResonantAbsorptionObservation(
        context=context,
        observation_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        pre_normalised=image,
    )
    bundle = build_rai_transmission(observation)
    np.testing.assert_array_equal(bundle.transmission, image.values)
    assert bundle.representation == "pre_normalised"
    assert any("outside [0, 1]" in item for item in bundle.warnings)

    wrong_quantity = ResonantAbsorptionObservation(
        context=context,
        observation_index=1,
        timestamp_s=1.0,
        pulse=pulse,
        pre_normalised=PreNormalisedCameraImage(
            np.ones((2, 2)),
            quantity_name="optical depth",
            normalisation_description="external processing",
        ),
    )
    with pytest.raises(ValueError, match="explicitly declare transmission"):
        build_rai_transmission(wrong_quantity)


def test_dual_port_binding_orders_h_then_v_and_flattens_roi(
    pulse: ProbePulse,
    dual_measurement: DualPortFaradayMeasurement,
) -> None:
    exposure = NonDestructiveExposure(
        exposure_index=3,
        timestamp_s=0.2,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"V": frame(200.0), "H": frame(100.0)},
    )
    bound = bind_non_destructive_exposure(
        exposure,
        dual_measurement,
        {"camera-1": calibration()},
    )
    assert tuple(bound.channel_electrons) == ("H", "V")
    np.testing.assert_array_equal(bound.channel_electrons["H"], 100.0)
    np.testing.assert_array_equal(bound.channel_electrons["V"], 200.0)
    pixel_count = dual_measurement.roi_pixel_count
    np.testing.assert_array_equal(bound.observed_vector_electrons[:pixel_count], 100.0)
    np.testing.assert_array_equal(bound.observed_vector_electrons[pixel_count:], 200.0)
    assert any("kappa_F" in item for item in bound.warnings)
    assert any("No TF, eGPE" in item for item in bound.assumptions)


def test_dark_field_binding_accepts_only_dark_field_channel(
    pulse: ProbePulse,
    dark_measurement: DarkFieldFaradayMeasurement,
) -> None:
    exposure = NonDestructiveExposure(
        exposure_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        readout="dark_field",
        raw_channels={"dark_field": frame(20.0, unit="adu")},
    )
    bound = bind_non_destructive_exposure(
        exposure,
        dark_measurement,
        {"camera-1": calibration(gain=0.5, bias=4.0)},
    )
    assert tuple(bound.channel_electrons) == ("dark_field",)
    np.testing.assert_allclose(bound.channel_electrons["dark_field"], 8.0)
    np.testing.assert_allclose(bound.observed_vector_electrons, 8.0)


def test_binding_rejects_readout_measurement_and_shape_mismatches(
    pulse: ProbePulse,
    dual_measurement: DualPortFaradayMeasurement,
    dark_measurement: DarkFieldFaradayMeasurement,
) -> None:
    dual_exposure = NonDestructiveExposure(
        exposure_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"H": frame(1.0), "V": frame(2.0)},
    )
    with pytest.raises(ValueError, match="incompatible"):
        bind_non_destructive_exposure(
            dual_exposure,
            dark_measurement,
            {"camera-1": calibration()},
        )

    wrong_shape = NonDestructiveExposure(
        exposure_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        readout="dual_port",
        raw_channels={
            "H": frame(1.0, shape=(3, 4)),
            "V": frame(2.0, shape=(3, 4)),
        },
    )
    with pytest.raises(ValueError, match="not measurement camera shape"):
        bind_non_destructive_exposure(
            wrong_shape,
            dual_measurement,
            {"camera-1": calibration()},
        )


def test_binding_requires_exact_camera_calibration_and_read_noise(
    pulse: ProbePulse,
    dual_measurement: DualPortFaradayMeasurement,
) -> None:
    exposure = NonDestructiveExposure(
        exposure_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"H": frame(1.0), "V": frame(2.0)},
    )
    with pytest.raises(ValueError, match="missing=.*camera-1"):
        bind_non_destructive_exposure(exposure, dual_measurement, {})
    with pytest.raises(ValueError, match="extra=.*camera-2"):
        bind_non_destructive_exposure(
            exposure,
            dual_measurement,
            {
                "camera-1": calibration(),
                "camera-2": calibration(camera_id="camera-2"),
            },
        )
    with pytest.raises(ValueError, match="read noise .* incompatible"):
        bind_non_destructive_exposure(
            exposure,
            dual_measurement,
            {"camera-1": calibration(read_noise=2.9)},
        )
    with pytest.raises(ValueError, match="mapping key"):
        bind_non_destructive_exposure(
            exposure,
            dual_measurement,
            {"camera-1": calibration(camera_id="camera-2")},
        )
