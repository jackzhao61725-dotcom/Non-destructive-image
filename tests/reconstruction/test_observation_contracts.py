"""Focused tests for morphology-free experimental input contracts."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, fields

import numpy as np
import pytest

from non_destructive_image.reconstruction.observations import (
    ExperimentContext,
    NonDestructiveExposure,
    NonDestructiveSequence,
    PreNormalisedCameraImage,
    ProbePulse,
    RawCameraFrame,
    ResonantAbsorptionObservation,
    UncertainQuantity,
)


@pytest.fixture
def context() -> ExperimentContext:
    return ExperimentContext(
        context_id="run-2026-07-21",
        condensate_quantities={
            "condensate_population": UncertainQuantity(2.5e4, 1.0e3, "1")
        },
        preparation_quantities={
            "bias_field": UncertainQuantity(1.2e-4, 2.0e-6, "T")
        },
        instrument_quantities={
            "magnification": UncertainQuantity(4.0, 0.1, "1"),
            "camera_pixel_pitch": UncertainQuantity(5.86e-6, 0.02e-6, "m"),
        },
        labels={"spin_preparation": "declared laboratory state"},
    )


@pytest.fixture
def pulse() -> ProbePulse:
    return ProbePulse(
        pulse_id="pulse-0",
        wavelength_m=401e-9,
        detuning_hz=-1.5e9,
        power_w=1.0e-3,
        duration_s=90e-6,
        polarisation_label="linear",
    )


def frame(
    value: float,
    *,
    shape: tuple[int, int] = (4, 5),
    unit: str = "electron",
    camera_id: str = "camera-1",
) -> RawCameraFrame:
    return RawCameraFrame(
        values=np.full(shape, value),
        unit=unit,  # type: ignore[arg-type]
        camera_id=camera_id,
    )


def exposure(
    index: int,
    timestamp_s: float,
    pulse: ProbePulse,
    *,
    h_value: float = 100.0,
) -> NonDestructiveExposure:
    return NonDestructiveExposure(
        exposure_index=index,
        timestamp_s=timestamp_s,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"H": frame(h_value), "V": frame(200.0)},
    )


def test_context_is_immutable_image_free_metadata_with_uncertainties(
    context: ExperimentContext,
) -> None:
    assert context.condensate_quantities["condensate_population"].unit == "1"
    assert not hasattr(context, "image")
    assert not hasattr(context, "reconstruct")
    with pytest.raises(TypeError):
        context.instrument_quantities["new"] = UncertainQuantity(  # type: ignore[index]
            1.0, 0.1, "1"
        )
    with pytest.raises(FrozenInstanceError):
        context.context_id = "changed"  # type: ignore[misc]


def test_quantity_requires_explicit_finite_nonnegative_uncertainty() -> None:
    with pytest.raises(ValueError, match="cannot be negative"):
        UncertainQuantity(1.0, -0.1, "m")
    with pytest.raises(ValueError, match="must be finite"):
        UncertainQuantity(1.0, float("nan"), "m")
    with pytest.raises(ValueError, match="cannot be empty"):
        UncertainQuantity(1.0, 0.1, "")


def test_raw_camera_arrays_are_validated_copied_and_read_only() -> None:
    source = np.arange(20.0).reshape(4, 5)
    recorded = RawCameraFrame(source, "adu", "camera-1")
    source[0, 0] = -999.0
    assert recorded.values[0, 0] == 0.0
    assert not recorded.values.flags.writeable
    with pytest.raises(ValueError, match="read-only"):
        recorded.values[0, 0] = 1.0
    with pytest.raises(ValueError, match="WRITEABLE"):
        recorded.values.setflags(write=True)
    with pytest.raises(ValueError, match="two-dimensional"):
        RawCameraFrame(np.ones(5), "adu", "camera-1")
    with pytest.raises(ValueError, match="finite"):
        RawCameraFrame(np.asarray([[np.nan]]), "adu", "camera-1")
    with pytest.raises(TypeError, match="real"):
        RawCameraFrame(np.asarray([[1.0 + 1.0j]]), "adu", "camera-1")


def test_rai_accepts_complete_raw_triplet(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    observation = ResonantAbsorptionObservation(
        context=context,
        observation_index=0,
        timestamp_s=1.25,
        pulse=pulse,
        atom_frame=frame(80.0),
        reference_frame=frame(100.0),
        dark_frame=frame(5.0),
    )
    assert observation.pre_normalised is None
    assert observation.atom_frame is not None
    assert observation.atom_frame.values.shape == (4, 5)


def test_rai_rejects_partial_or_mixed_representations(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    normalised = PreNormalisedCameraImage(
        np.full((4, 5), 0.8),
        quantity_name="transmission",
        normalisation_description="(atom-dark)/(reference-dark)",
    )
    with pytest.raises(ValueError, match="requires atom, reference and dark"):
        ResonantAbsorptionObservation(
            context=context,
            observation_index=0,
            timestamp_s=0.0,
            pulse=pulse,
            atom_frame=frame(80.0),
        )
    with pytest.raises(ValueError, match="exactly one"):
        ResonantAbsorptionObservation(
            context=context,
            observation_index=0,
            timestamp_s=0.0,
            pulse=pulse,
            atom_frame=frame(80.0),
            reference_frame=frame(100.0),
            dark_frame=frame(5.0),
            pre_normalised=normalised,
        )
    accepted = ResonantAbsorptionObservation(
        context=context,
        observation_index=0,
        timestamp_s=0.0,
        pulse=pulse,
        pre_normalised=normalised,
    )
    assert accepted.pre_normalised is normalised


@pytest.mark.parametrize(
    "reference,dark,match",
    [
        (frame(100.0, shape=(3, 5)), frame(5.0), "share a shape"),
        (frame(100.0, unit="adu"), frame(5.0), "share a unit"),
        (frame(100.0, camera_id="camera-2"), frame(5.0), "share a camera_id"),
    ],
)
def test_rai_raw_triplet_requires_one_camera_contract(
    context: ExperimentContext,
    pulse: ProbePulse,
    reference: RawCameraFrame,
    dark: RawCameraFrame,
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        ResonantAbsorptionObservation(
            context=context,
            observation_index=0,
            timestamp_s=0.0,
            pulse=pulse,
            atom_frame=frame(80.0),
            reference_frame=reference,
            dark_frame=dark,
        )


def test_dual_port_exposure_records_simultaneous_h_and_v(
    pulse: ProbePulse,
) -> None:
    recorded = exposure(0, 0.0, pulse)
    assert tuple(recorded.raw_channels) == ("H", "V")
    with pytest.raises(ValueError, match="exactly simultaneous H and V"):
        NonDestructiveExposure(
            exposure_index=0,
            timestamp_s=0.0,
            pulse=pulse,
            readout="dual_port",
            raw_channels={"H": frame(100.0)},
        )
    with pytest.raises(ValueError, match="share a shape"):
        NonDestructiveExposure(
            exposure_index=0,
            timestamp_s=0.0,
            pulse=pulse,
            readout="dual_port",
            raw_channels={"H": frame(100.0), "V": frame(200.0, shape=(3, 5))},
        )


def test_sequence_validates_index_time_and_pulse_order(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    first = exposure(4, 0.0, pulse)
    second = exposure(7, 1.0e-3, pulse)
    sequence = NonDestructiveSequence("sequence-1", context, (first, second))
    assert [item.exposure_index for item in sequence.exposures] == [4, 7]

    with pytest.raises(ValueError, match="indices must increase"):
        NonDestructiveSequence("bad-index", context, (first, exposure(3, 1.0e-3, pulse)))
    with pytest.raises(ValueError, match="timestamps must increase"):
        NonDestructiveSequence("bad-time", context, (first, exposure(7, 0.0, pulse)))
    with pytest.raises(ValueError, match="cannot overlap"):
        NonDestructiveSequence("overlap", context, (first, exposure(7, 40e-6, pulse)))


def test_sequence_rejects_camera_contract_changes(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    first = exposure(0, 0.0, pulse)
    changed = NonDestructiveExposure(
        exposure_index=1,
        timestamp_s=1.0e-3,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"H": frame(100.0, camera_id="camera-2"), "V": frame(200.0)},
    )
    with pytest.raises(ValueError, match="camera shape, unit and camera_id"):
        NonDestructiveSequence("changed-camera", context, (first, changed))


def test_sequence_channel_mapping_order_is_not_physical_order(
    context: ExperimentContext,
    pulse: ProbePulse,
) -> None:
    first = exposure(0, 0.0, pulse)
    reversed_mapping = NonDestructiveExposure(
        exposure_index=1,
        timestamp_s=1.0e-3,
        pulse=pulse,
        readout="dual_port",
        raw_channels={"V": frame(200.0), "H": frame(100.0)},
    )
    sequence = NonDestructiveSequence(
        "mapping-order",
        context,
        (first, reversed_mapping),
    )
    assert len(sequence.exposures) == 2


def test_contract_field_names_contain_no_morphology_model_assumptions() -> None:
    contract_types = (
        ExperimentContext,
        ResonantAbsorptionObservation,
        NonDestructiveExposure,
        NonDestructiveSequence,
    )
    names = {
        definition.name.lower()
        for contract in contract_types
        for definition in fields(contract)
    }
    forbidden = ("thomas", "fermi", "egpe", "morphology", "density_model")
    assert not any(token in name for name in names for token in forbidden)
