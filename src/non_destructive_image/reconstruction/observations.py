"""Immutable input contracts for experimental reconstruction data.

The classes in this module describe what was prepared, what probe pulse was
applied and what the camera recorded.  They deliberately contain no object
model, morphology label or reconstruction method.  An :class:`ExperimentContext`
therefore cannot reconstruct anything by itself; it must be paired with one or
more observations and an independently selected inverse model.

Numeric physical quantities use SI units named by their field suffixes.  Camera
arrays are either raw analogue-to-digital units (ADU), calibrated electron
counts, or an explicitly documented dimensionless pre-normalised image.
Timestamps are seconds on the experiment clock and denote pulse start times.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


RawCameraUnit = Literal["adu", "electron"]


def _nonempty_text(value: str, label: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} cannot be empty")
    return text


def _finite_scalar(value: float, label: str) -> float:
    scalar = float(value)
    if not np.isfinite(scalar):
        raise ValueError(f"{label} must be finite")
    return scalar


def _frozen_string_mapping(
    values: Mapping[str, str],
    label: str,
) -> Mapping[str, str]:
    copied: dict[str, str] = {}
    for key, value in values.items():
        name = _nonempty_text(key, f"{label} key")
        if name in copied:
            raise ValueError(f"{label} contains a duplicate key: {name!r}")
        copied[name] = _nonempty_text(value, f"{label}[{name!r}]")
    return MappingProxyType(copied)


def _frozen_quantity_mapping(
    values: Mapping[str, "UncertainQuantity"],
    label: str,
) -> Mapping[str, "UncertainQuantity"]:
    copied: dict[str, UncertainQuantity] = {}
    for key, value in values.items():
        name = _nonempty_text(key, f"{label} key")
        if not isinstance(value, UncertainQuantity):
            raise TypeError(f"{label}[{name!r}] must be an UncertainQuantity")
        if name in copied:
            raise ValueError(f"{label} contains a duplicate key: {name!r}")
        copied[name] = value
    if not copied:
        raise ValueError(f"{label} cannot be empty")
    return MappingProxyType(copied)


def _frozen_camera_array(values: ArrayLike, label: str) -> NDArray[np.floating]:
    array = np.asarray(values)
    if array.ndim != 2 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty two-dimensional array")
    if np.issubdtype(array.dtype, np.bool_) or not np.issubdtype(
        array.dtype, np.number
    ):
        raise TypeError(f"{label} must contain real numeric camera values")
    if np.iscomplexobj(array):
        raise TypeError(f"{label} must contain real camera values")
    copied = np.asarray(array, dtype=float).copy()
    if np.any(~np.isfinite(copied)):
        raise ValueError(f"{label} must contain only finite camera values")
    frozen = np.frombuffer(copied.tobytes(order="C"), dtype=copied.dtype).reshape(
        copied.shape
    )
    return frozen


@dataclass(frozen=True)
class UncertainQuantity:
    """One scalar value and its one-standard-deviation uncertainty.

    ``value_si`` and ``standard_uncertainty_si`` are expressed in the SI unit
    named by ``unit``.  Use ``"1"`` for dimensionless quantities.  A zero
    uncertainty is permitted for an exact declared setpoint, but missing or
    unknown uncertainties must not be silently represented by NaN.
    """

    value_si: float
    standard_uncertainty_si: float
    unit: str
    source: str = "declared"

    def __post_init__(self) -> None:
        value = _finite_scalar(self.value_si, "quantity value")
        uncertainty = _finite_scalar(
            self.standard_uncertainty_si, "quantity standard uncertainty"
        )
        if uncertainty < 0.0:
            raise ValueError("quantity standard uncertainty cannot be negative")
        object.__setattr__(self, "value_si", value)
        object.__setattr__(self, "standard_uncertainty_si", uncertainty)
        object.__setattr__(self, "unit", _nonempty_text(self.unit, "quantity unit"))
        object.__setattr__(self, "source", _nonempty_text(self.source, "quantity source"))


@dataclass(frozen=True)
class ExperimentContext:
    """Shared preparation and instrument metadata, without image data.

    The three quantity mappings are intentionally generic.  Their keys describe
    measured or declared scalar metadata and their values always carry an
    uncertainty and unit.  Optional labels hold categorical information such as
    a preparation state or camera serial number.  No TF, eGPE or other
    morphology assumption is encoded here.
    """

    context_id: str
    condensate_quantities: Mapping[str, UncertainQuantity]
    preparation_quantities: Mapping[str, UncertainQuantity]
    instrument_quantities: Mapping[str, UncertainQuantity]
    labels: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "context_id", _nonempty_text(self.context_id, "context_id"))
        object.__setattr__(
            self,
            "condensate_quantities",
            _frozen_quantity_mapping(self.condensate_quantities, "condensate_quantities"),
        )
        object.__setattr__(
            self,
            "preparation_quantities",
            _frozen_quantity_mapping(self.preparation_quantities, "preparation_quantities"),
        )
        object.__setattr__(
            self,
            "instrument_quantities",
            _frozen_quantity_mapping(self.instrument_quantities, "instrument_quantities"),
        )
        object.__setattr__(self, "labels", _frozen_string_mapping(self.labels, "labels"))


@dataclass(frozen=True)
class ProbePulse:
    """Probe-pulse setpoints in SI units.

    ``detuning_hz`` is signed and may be zero for resonant imaging.  The power
    is the delivered optical power associated with the observation and the
    duration is the pulse duration.
    """

    pulse_id: str
    wavelength_m: float
    detuning_hz: float
    power_w: float
    duration_s: float
    polarisation_label: str

    def __post_init__(self) -> None:
        wavelength = _finite_scalar(self.wavelength_m, "pulse wavelength_m")
        detuning = _finite_scalar(self.detuning_hz, "pulse detuning_hz")
        power = _finite_scalar(self.power_w, "pulse power_w")
        duration = _finite_scalar(self.duration_s, "pulse duration_s")
        if wavelength <= 0.0:
            raise ValueError("pulse wavelength_m must be positive")
        if power <= 0.0:
            raise ValueError("pulse power_w must be positive")
        if duration <= 0.0:
            raise ValueError("pulse duration_s must be positive")
        object.__setattr__(self, "pulse_id", _nonempty_text(self.pulse_id, "pulse_id"))
        object.__setattr__(self, "wavelength_m", wavelength)
        object.__setattr__(self, "detuning_hz", detuning)
        object.__setattr__(self, "power_w", power)
        object.__setattr__(self, "duration_s", duration)
        object.__setattr__(
            self,
            "polarisation_label",
            _nonempty_text(self.polarisation_label, "polarisation_label"),
        )


@dataclass(frozen=True)
class RawCameraFrame:
    """One raw camera plane in ADU or calibrated electron counts."""

    values: NDArray[np.floating] = field(repr=False, compare=False)
    unit: RawCameraUnit
    camera_id: str

    def __post_init__(self) -> None:
        if self.unit not in {"adu", "electron"}:
            raise ValueError("raw camera unit must be 'adu' or 'electron'")
        object.__setattr__(self, "values", _frozen_camera_array(self.values, "camera frame"))
        object.__setattr__(self, "camera_id", _nonempty_text(self.camera_id, "camera_id"))


@dataclass(frozen=True)
class PreNormalisedCameraImage:
    """Explicit dimensionless equivalent of an RAI atom/reference/dark triplet."""

    values: NDArray[np.floating] = field(repr=False, compare=False)
    quantity_name: str
    normalisation_description: str
    unit: Literal["1"] = "1"

    def __post_init__(self) -> None:
        if self.unit != "1":
            raise ValueError("a pre-normalised camera image must be dimensionless")
        object.__setattr__(
            self,
            "values",
            _frozen_camera_array(self.values, "pre-normalised camera image"),
        )
        object.__setattr__(
            self,
            "quantity_name",
            _nonempty_text(self.quantity_name, "pre-normalised quantity_name"),
        )
        object.__setattr__(
            self,
            "normalisation_description",
            _nonempty_text(
                self.normalisation_description,
                "pre-normalised normalisation_description",
            ),
        )


@dataclass(frozen=True)
class ResonantAbsorptionObservation:
    """One RAI observation with raw triplet or explicit pre-normalised data.

    Exactly one representation is accepted: either complete atom, reference and
    dark camera frames, or one :class:`PreNormalisedCameraImage`.  Partial raw
    triplets and mixtures of the two representations are rejected.
    """

    context: ExperimentContext
    observation_index: int
    timestamp_s: float
    pulse: ProbePulse
    atom_frame: RawCameraFrame | None = None
    reference_frame: RawCameraFrame | None = None
    dark_frame: RawCameraFrame | None = None
    pre_normalised: PreNormalisedCameraImage | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.context, ExperimentContext):
            raise TypeError("context must be an ExperimentContext")
        if not isinstance(self.observation_index, int) or isinstance(
            self.observation_index, bool
        ):
            raise TypeError("observation_index must be an integer")
        if self.observation_index < 0:
            raise ValueError("observation_index cannot be negative")
        timestamp = _finite_scalar(self.timestamp_s, "RAI timestamp_s")
        if timestamp < 0.0:
            raise ValueError("RAI timestamp_s cannot be negative")
        if not isinstance(self.pulse, ProbePulse):
            raise TypeError("pulse must be a ProbePulse")
        object.__setattr__(self, "timestamp_s", timestamp)

        raw = (self.atom_frame, self.reference_frame, self.dark_frame)
        raw_complete = all(frame is not None for frame in raw)
        raw_present = any(frame is not None for frame in raw)
        normalised_present = self.pre_normalised is not None
        if raw_present and not raw_complete:
            raise ValueError("RAI raw input requires atom, reference and dark frames")
        if raw_complete == normalised_present:
            raise ValueError(
                "RAI input requires exactly one of a raw triplet or a pre-normalised image"
            )
        if raw_complete:
            frames = tuple(frame for frame in raw if frame is not None)
            if not all(isinstance(frame, RawCameraFrame) for frame in frames):
                raise TypeError("RAI raw inputs must be RawCameraFrame instances")
            shapes = {frame.values.shape for frame in frames}
            units = {frame.unit for frame in frames}
            cameras = {frame.camera_id for frame in frames}
            if len(shapes) != 1:
                raise ValueError("RAI atom, reference and dark frames must share a shape")
            if len(units) != 1:
                raise ValueError("RAI atom, reference and dark frames must share a unit")
            if len(cameras) != 1:
                raise ValueError("RAI atom, reference and dark frames must share a camera_id")
        elif not isinstance(self.pre_normalised, PreNormalisedCameraImage):
            raise TypeError("pre_normalised must be a PreNormalisedCameraImage")


@dataclass(frozen=True)
class NonDestructiveExposure:
    """One ordered non-destructive exposure and its simultaneous raw channels.

    A ``readout`` value of ``"dual_port"`` requires exactly ``H`` and ``V``.
    Because both frames belong to this single exposure and timestamp, the
    contract records them as a simultaneous pair rather than two sequential
    observations.
    """

    exposure_index: int
    timestamp_s: float
    pulse: ProbePulse
    readout: str
    raw_channels: Mapping[str, RawCameraFrame]

    def __post_init__(self) -> None:
        if not isinstance(self.exposure_index, int) or isinstance(
            self.exposure_index, bool
        ):
            raise TypeError("exposure_index must be an integer")
        if self.exposure_index < 0:
            raise ValueError("exposure_index cannot be negative")
        timestamp = _finite_scalar(self.timestamp_s, "exposure timestamp_s")
        if timestamp < 0.0:
            raise ValueError("exposure timestamp_s cannot be negative")
        if not isinstance(self.pulse, ProbePulse):
            raise TypeError("pulse must be a ProbePulse")
        readout = _nonempty_text(self.readout, "readout")
        channels: dict[str, RawCameraFrame] = {}
        for key, frame in self.raw_channels.items():
            name = _nonempty_text(key, "raw channel name")
            if not isinstance(frame, RawCameraFrame):
                raise TypeError(f"raw channel {name!r} must be a RawCameraFrame")
            if name in channels:
                raise ValueError(f"duplicate raw channel name: {name!r}")
            channels[name] = frame
        if not channels:
            raise ValueError("an exposure must contain at least one raw camera channel")
        if readout == "dual_port":
            if set(channels) != {"H", "V"}:
                raise ValueError(
                    "dual_port exposure requires exactly simultaneous H and V channels"
                )
            h_frame, v_frame = channels["H"], channels["V"]
            if h_frame.values.shape != v_frame.values.shape:
                raise ValueError("dual_port H and V channels must share a shape")
            if h_frame.unit != v_frame.unit:
                raise ValueError("dual_port H and V channels must share a unit")
        object.__setattr__(self, "timestamp_s", timestamp)
        object.__setattr__(self, "readout", readout)
        object.__setattr__(self, "raw_channels", MappingProxyType(channels))


@dataclass(frozen=True)
class NonDestructiveSequence:
    """Ordered, non-overlapping exposures acquired under one context.

    Indices and pulse-start timestamps must increase strictly.  The next pulse
    may start when the preceding pulse ends, but not before it.  Channel names,
    shapes, units and camera assignments remain fixed across the sequence so a
    sequence cannot silently mix incompatible camera contracts.
    """

    sequence_id: str
    context: ExperimentContext
    exposures: tuple[NonDestructiveExposure, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.context, ExperimentContext):
            raise TypeError("context must be an ExperimentContext")
        sequence_id = _nonempty_text(self.sequence_id, "sequence_id")
        exposures = tuple(self.exposures)
        if not exposures:
            raise ValueError("a non-destructive sequence cannot be empty")
        if not all(isinstance(exposure, NonDestructiveExposure) for exposure in exposures):
            raise TypeError("all sequence entries must be NonDestructiveExposure instances")

        first = exposures[0]
        expected_readout = first.readout
        expected_channels = tuple(sorted(first.raw_channels))
        expected_contract = {
            name: (
                first.raw_channels[name].values.shape,
                first.raw_channels[name].unit,
                first.raw_channels[name].camera_id,
            )
            for name in expected_channels
        }
        for previous, current in zip(exposures[:-1], exposures[1:], strict=True):
            if current.exposure_index <= previous.exposure_index:
                raise ValueError("sequence exposure indices must increase strictly")
            if current.timestamp_s <= previous.timestamp_s:
                raise ValueError("sequence timestamps must increase strictly")
            previous_end = previous.timestamp_s + previous.pulse.duration_s
            if current.timestamp_s < previous_end:
                raise ValueError("sequence probe pulses cannot overlap")
        for exposure in exposures:
            if exposure.readout != expected_readout:
                raise ValueError("sequence readout must remain constant")
            if set(exposure.raw_channels) != set(expected_channels):
                raise ValueError("sequence raw channel names must remain constant")
            actual_contract = {
                name: (
                    exposure.raw_channels[name].values.shape,
                    exposure.raw_channels[name].unit,
                    exposure.raw_channels[name].camera_id,
                )
                for name in expected_channels
            }
            if actual_contract != expected_contract:
                raise ValueError("sequence camera shape, unit and camera_id must remain constant")

        object.__setattr__(self, "sequence_id", sequence_id)
        object.__setattr__(self, "exposures", exposures)
