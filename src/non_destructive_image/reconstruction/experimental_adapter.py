"""Strict adapters from experimental camera contracts to inverse inputs.

This module performs only camera-level transformations whose meaning is
declared by the input contracts in :mod:`.observations`.  It does not convert
absorption transmission to density, choose an object morphology, or infer the
Faraday response coefficient.  Those operations belong to separately stated
physical and inverse models.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .measurements import DarkFieldFaradayMeasurement, DualPortFaradayMeasurement
from .observations import (
    NonDestructiveExposure,
    RawCameraFrame,
    ResonantAbsorptionObservation,
)


FaradayMeasurement = DualPortFaradayMeasurement | DarkFieldFaradayMeasurement


def _nonempty_text(value: str, label: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} cannot be empty")
    return text


def _finite_nonnegative(value: float, label: str) -> float:
    scalar = float(value)
    if not np.isfinite(scalar):
        raise ValueError(f"{label} must be finite")
    if scalar < 0.0:
        raise ValueError(f"{label} cannot be negative")
    return scalar


def _frozen_float_array(values: ArrayLike, label: str) -> NDArray[np.floating]:
    array = np.asarray(values)
    if array.ndim != 2 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty two-dimensional array")
    if np.issubdtype(array.dtype, np.bool_) or not np.issubdtype(array.dtype, np.number):
        raise TypeError(f"{label} must contain real numeric values")
    if np.iscomplexobj(array):
        raise TypeError(f"{label} must contain real values")
    copied = np.asarray(array, dtype=float).copy()
    if np.any(~np.isfinite(copied)):
        raise ValueError(f"{label} must contain only finite values")
    return np.frombuffer(copied.tobytes(order="C"), dtype=float).reshape(copied.shape)


def _frozen_bool_array(values: ArrayLike, label: str) -> NDArray[np.bool_]:
    array = np.asarray(values)
    if array.ndim != 2 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty two-dimensional array")
    copied = np.asarray(array, dtype=bool).copy()
    return np.frombuffer(copied.tobytes(order="C"), dtype=bool).reshape(copied.shape)


def _frozen_vector(values: ArrayLike, label: str) -> NDArray[np.floating]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty one-dimensional array")
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{label} must contain only finite values")
    copied = array.copy()
    return np.frombuffer(copied.tobytes(order="C"), dtype=float)


@dataclass(frozen=True)
class CameraCalibration:
    """Camera conversion and read-noise contract for one camera.

    ``bias_adu`` may be a scalar or a complete two-dimensional bias map.  ADU
    frames are converted as ``(adu - bias_adu) * electrons_per_adu``.  Frames
    already labelled ``"electron"`` are treated as calibrated electron counts:
    neither the ADU gain nor bias is applied to them.  The read-noise value is
    retained so that the adapter can check compatibility with the measurement
    operator's Poisson--Gaussian likelihood.
    """

    camera_id: str
    electrons_per_adu: float
    read_noise_electrons: float
    bias_adu: float | NDArray[np.floating] = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        camera_id = _nonempty_text(self.camera_id, "camera_id")
        gain = float(self.electrons_per_adu)
        if not np.isfinite(gain) or gain <= 0.0:
            raise ValueError("electrons_per_adu must be finite and positive")
        read_noise = _finite_nonnegative(
            self.read_noise_electrons,
            "read_noise_electrons",
        )
        bias = np.asarray(self.bias_adu)
        if bias.ndim == 0:
            scalar_bias = float(bias)
            if not np.isfinite(scalar_bias):
                raise ValueError("scalar bias_adu must be finite")
            frozen_bias: float | NDArray[np.floating] = scalar_bias
        else:
            frozen_bias = _frozen_float_array(bias, "bias_adu map")
        object.__setattr__(self, "camera_id", camera_id)
        object.__setattr__(self, "electrons_per_adu", gain)
        object.__setattr__(self, "read_noise_electrons", read_noise)
        object.__setattr__(self, "bias_adu", frozen_bias)

    def bias_for_shape(self, shape: tuple[int, int]) -> float | NDArray[np.floating]:
        """Return the scalar/map bias after checking a frame shape."""

        if isinstance(self.bias_adu, np.ndarray) and self.bias_adu.shape != shape:
            raise ValueError(
                f"camera {self.camera_id!r} bias map has shape {self.bias_adu.shape}, "
                f"not frame shape {shape}"
            )
        return self.bias_adu


@dataclass(frozen=True)
class CalibratedCameraFrame:
    """One immutable camera plane expressed in detected electrons."""

    values_electrons: NDArray[np.floating] = field(repr=False, compare=False)
    camera_id: str
    source_unit: Literal["adu", "electron"]
    assumptions: tuple[str, ...]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "values_electrons",
            _frozen_float_array(self.values_electrons, "calibrated electron frame"),
        )
        object.__setattr__(self, "camera_id", _nonempty_text(self.camera_id, "camera_id"))
        if self.source_unit not in {"adu", "electron"}:
            raise ValueError("source_unit must be 'adu' or 'electron'")
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        object.__setattr__(self, "warnings", tuple(self.warnings))


@dataclass(frozen=True)
class TransmissionBundle:
    """Camera-level RAI transmission without a density interpretation."""

    transmission: NDArray[np.floating] = field(repr=False, compare=False)
    representation: Literal["raw_triplet", "pre_normalised"]
    valid_mask: NDArray[np.bool_] = field(repr=False, compare=False)
    atom_electrons: NDArray[np.floating] | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    reference_electrons: NDArray[np.floating] | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    dark_electrons: NDArray[np.floating] | None = field(
        default=None,
        repr=False,
        compare=False,
    )
    assumptions: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        transmission = _frozen_float_array(self.transmission, "transmission")
        valid_mask = _frozen_bool_array(self.valid_mask, "transmission valid_mask")
        if transmission.shape != valid_mask.shape:
            raise ValueError("transmission and valid_mask must share a shape")
        object.__setattr__(self, "transmission", transmission)
        object.__setattr__(self, "valid_mask", valid_mask)
        for name in ("atom_electrons", "reference_electrons", "dark_electrons"):
            values = getattr(self, name)
            if values is not None:
                frozen = _frozen_float_array(values, name)
                if frozen.shape != transmission.shape:
                    raise ValueError(f"{name} must share the transmission shape")
                object.__setattr__(self, name, frozen)
        if self.representation == "raw_triplet":
            if any(
                getattr(self, name) is None
                for name in ("atom_electrons", "reference_electrons", "dark_electrons")
            ):
                raise ValueError("raw_triplet transmission must retain all three electron frames")
        elif self.representation == "pre_normalised":
            if any(
                getattr(self, name) is not None
                for name in ("atom_electrons", "reference_electrons", "dark_electrons")
            ):
                raise ValueError("pre_normalised transmission cannot contain raw electron frames")
        else:
            raise ValueError("unknown transmission representation")
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        object.__setattr__(self, "warnings", tuple(self.warnings))


@dataclass(frozen=True)
class BoundNonDestructiveExposure:
    """Electron-count channels bound to one compatible Faraday operator."""

    exposure_index: int
    timestamp_s: float
    readout: Literal["dual_port", "dark_field"]
    channel_electrons: Mapping[str, NDArray[np.floating]] = field(
        repr=False,
        compare=False,
    )
    observed_vector_electrons: NDArray[np.floating] = field(
        repr=False,
        compare=False,
    )
    assumptions: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        frozen_channels: dict[str, NDArray[np.floating]] = {}
        for name, values in self.channel_electrons.items():
            channel_name = _nonempty_text(name, "channel name")
            frozen_channels[channel_name] = _frozen_float_array(
                values,
                f"{channel_name} electron channel",
            )
        expected = ("H", "V") if self.readout == "dual_port" else ("dark_field",)
        if tuple(frozen_channels) != expected:
            raise ValueError(f"{self.readout} channel order must be {expected}")
        object.__setattr__(self, "channel_electrons", MappingProxyType(frozen_channels))
        object.__setattr__(
            self,
            "observed_vector_electrons",
            _frozen_vector(self.observed_vector_electrons, "observed electron vector"),
        )
        object.__setattr__(self, "assumptions", tuple(self.assumptions))
        object.__setattr__(self, "warnings", tuple(self.warnings))


def calibrate_camera_frame(
    frame: RawCameraFrame,
    calibration: CameraCalibration | None,
) -> CalibratedCameraFrame:
    """Convert one frame to electrons without modifying electron-labelled data."""

    if not isinstance(frame, RawCameraFrame):
        raise TypeError("frame must be a RawCameraFrame")
    if frame.unit == "adu":
        if calibration is None:
            raise ValueError("ADU input requires an explicit CameraCalibration")
        if frame.camera_id != calibration.camera_id:
            raise ValueError(
                f"frame camera_id {frame.camera_id!r} does not match calibration "
                f"{calibration.camera_id!r}"
            )
        bias = calibration.bias_for_shape(frame.values.shape)
        values = (frame.values - bias) * calibration.electrons_per_adu
        assumptions = (
            "The declared bias and electrons-per-ADU gain apply to this complete frame.",
        )
        warnings: tuple[str, ...] = ()
    else:
        if calibration is not None:
            if frame.camera_id != calibration.camera_id:
                raise ValueError(
                    f"frame camera_id {frame.camera_id!r} does not match calibration "
                    f"{calibration.camera_id!r}"
                )
        values = frame.values
        assumptions = (
            "Electron-labelled input is already bias- and gain-calibrated; no ADU conversion is applied.",
        )
        warnings = (
            "The adapter cannot independently verify the upstream electron calibration.",
        )
    return CalibratedCameraFrame(
        values_electrons=values,
        camera_id=frame.camera_id,
        source_unit=frame.unit,
        assumptions=assumptions,
        warnings=warnings,
    )


def build_rai_transmission(
    observation: ResonantAbsorptionObservation,
    calibration: CameraCalibration | None = None,
) -> TransmissionBundle:
    """Build the declared RAI transmission and stop before density inversion."""

    if not isinstance(observation, ResonantAbsorptionObservation):
        raise TypeError("observation must be a ResonantAbsorptionObservation")
    if observation.pre_normalised is not None:
        quantity = "".join(
            character
            for character in observation.pre_normalised.quantity_name.lower()
            if character.isalnum()
        )
        if quantity not in {
            "transmission",
            "normalisedtransmission",
            "normalizedtransmission",
        }:
            raise ValueError(
                "pre-normalised RAI input must explicitly declare transmission as its quantity"
            )
        values = observation.pre_normalised.values
        warnings: list[str] = [
            "The adapter accepts the supplied normalisation and cannot reconstruct its raw-frame covariance."
        ]
        if np.any((values < 0.0) | (values > 1.0)):
            warnings.append(
                "Transmission contains values outside [0, 1]; they are retained rather than clipped."
            )
        return TransmissionBundle(
            transmission=values,
            representation="pre_normalised",
            valid_mask=np.ones(values.shape, dtype=bool),
            assumptions=(
                observation.pre_normalised.normalisation_description,
                "The dimensionless input is interpreted only as camera-level transmission.",
                "No Beer--Lambert, saturation, TF or eGPE density model is applied.",
            ),
            warnings=tuple(warnings),
        )

    atom = calibrate_camera_frame(observation.atom_frame, calibration)  # type: ignore[arg-type]
    reference = calibrate_camera_frame(observation.reference_frame, calibration)  # type: ignore[arg-type]
    dark = calibrate_camera_frame(observation.dark_frame, calibration)  # type: ignore[arg-type]
    denominator = reference.values_electrons - dark.values_electrons
    if np.any(denominator <= 0.0):
        count = int(np.count_nonzero(denominator <= 0.0))
        raise ValueError(
            f"RAI reference-minus-dark must be positive in every pixel; {count} pixels fail"
        )
    transmission = (atom.values_electrons - dark.values_electrons) / denominator
    warnings = list(atom.warnings + reference.warnings + dark.warnings)
    if np.any((transmission < 0.0) | (transmission > 1.0)):
        warnings.append(
            "Transmission contains values outside [0, 1]; they are retained rather than clipped."
        )
    return TransmissionBundle(
        transmission=transmission,
        representation="raw_triplet",
        valid_mask=np.ones(transmission.shape, dtype=bool),
        atom_electrons=atom.values_electrons,
        reference_electrons=reference.values_electrons,
        dark_electrons=dark.values_electrons,
        assumptions=(
            "Transmission is (atom-dark)/(reference-dark) after the declared camera conversion.",
            "The three frames are registered pixel-for-pixel and share one camera response.",
            "No Beer--Lambert, saturation, TF or eGPE density model is applied.",
        ),
        warnings=tuple(dict.fromkeys(warnings)),
    )


def _calibrations_by_camera(
    calibrations: Mapping[str, CameraCalibration],
    required_camera_ids: set[str],
) -> Mapping[str, CameraCalibration]:
    copied: dict[str, CameraCalibration] = {}
    for key, calibration in calibrations.items():
        camera_id = _nonempty_text(key, "camera-calibration key")
        if not isinstance(calibration, CameraCalibration):
            raise TypeError(f"calibration {camera_id!r} must be a CameraCalibration")
        if camera_id != calibration.camera_id:
            raise ValueError(
                f"calibration mapping key {camera_id!r} does not match "
                f"CameraCalibration.camera_id {calibration.camera_id!r}"
            )
        copied[camera_id] = calibration
    if set(copied) != required_camera_ids:
        missing = sorted(required_camera_ids - set(copied))
        extra = sorted(set(copied) - required_camera_ids)
        raise ValueError(
            f"camera calibrations must match exposure cameras exactly; missing={missing}, extra={extra}"
        )
    return MappingProxyType(copied)


def bind_non_destructive_exposure(
    exposure: NonDestructiveExposure,
    measurement: FaradayMeasurement,
    calibrations: Mapping[str, CameraCalibration],
) -> BoundNonDestructiveExposure:
    """Validate, calibrate and order one exposure for a Faraday likelihood."""

    if not isinstance(exposure, NonDestructiveExposure):
        raise TypeError("exposure must be a NonDestructiveExposure")
    if isinstance(measurement, DualPortFaradayMeasurement):
        expected_readout: Literal["dual_port", "dark_field"] = "dual_port"
        expected_channels = ("H", "V")
    elif isinstance(measurement, DarkFieldFaradayMeasurement):
        expected_readout = "dark_field"
        expected_channels = ("dark_field",)
    else:
        raise TypeError(
            "measurement must be a DualPortFaradayMeasurement or "
            "DarkFieldFaradayMeasurement"
        )
    if exposure.readout != expected_readout:
        raise ValueError(
            f"exposure readout {exposure.readout!r} is incompatible with "
            f"{type(measurement).__name__}"
        )
    if set(exposure.raw_channels) != set(expected_channels):
        raise ValueError(
            f"{expected_readout} exposure requires exactly the channels {expected_channels}"
        )
    required_camera_ids = {
        exposure.raw_channels[name].camera_id for name in expected_channels
    }
    camera_calibrations = _calibrations_by_camera(calibrations, required_camera_ids)
    calibrated: dict[str, NDArray[np.floating]] = {}
    assumptions: list[str] = []
    warnings: list[str] = []
    for name in expected_channels:
        raw = exposure.raw_channels[name]
        if raw.values.shape != measurement.camera_shape:
            raise ValueError(
                f"channel {name!r} has shape {raw.values.shape}, not measurement "
                f"camera shape {measurement.camera_shape}"
            )
        calibration = camera_calibrations[raw.camera_id]
        if not np.isclose(
            calibration.read_noise_electrons,
            measurement.read_noise_electrons,
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                f"camera {raw.camera_id!r} read noise "
                f"{calibration.read_noise_electrons} e- is incompatible with measurement "
                f"value {measurement.read_noise_electrons} e-"
            )
        converted = calibrate_camera_frame(raw, calibration)
        calibrated[name] = converted.values_electrons
        assumptions.extend(converted.assumptions)
        warnings.extend(converted.warnings)
    observed_vector = measurement.flatten_observed(
        *(calibrated[name] for name in expected_channels)
    )
    assumptions.extend(
        [
            "Channel ordering follows the measurement operator rather than mapping insertion order.",
            "The measurement operator's fixed aperture, detector scale and Faraday response are retained.",
            "No TF, eGPE or other condensate morphology is introduced by this adapter.",
        ]
    )
    warnings.extend(
        [
            "The adapter does not infer or calibrate kappa_F from the camera data.",
            "Absolute density remains conditional on the independently supplied Faraday response contract.",
        ]
    )
    return BoundNonDestructiveExposure(
        exposure_index=exposure.exposure_index,
        timestamp_s=exposure.timestamp_s,
        readout=expected_readout,
        channel_electrons=calibrated,
        observed_vector_electrons=observed_vector,
        assumptions=tuple(dict.fromkeys(assumptions)),
        warnings=tuple(dict.fromkeys(warnings)),
    )
