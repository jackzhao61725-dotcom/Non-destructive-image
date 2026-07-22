"""Analytic column-density truth inputs for reconstruction stress tests.

The providers in this module generate non-negative two-dimensional column
density maps.  They are deliberately independent of the fitting basis and the
measurement operator, so they can be used to expose model mismatch without
coupling a truth generator to the inverse model.

These maps are controlled analytic test patterns.  They are not solutions of a
Gross--Pitaevskii equation, predictions for dipolar erbium, or classifiers for
condensate phases.  Latent-state metadata may be attached to identical density
maps to make density-only impossibility tests explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .synthetic_morphologies import SyntheticMorphology


MetadataValue: TypeAlias = str | int | float | bool | tuple[float, ...]


@dataclass(frozen=True)
class AnalyticTruthInput:
    """One analytic truth map plus generation and latent-state metadata.

    ``observational_equivalence_group`` labels inputs that intentionally share
    the same column density while differing in metadata that a density-only
    measurement cannot observe.  It records an impossibility test; it does not
    imply that the optical model has measured the latent label.
    """

    morphology: SyntheticMorphology
    provider: str
    parameter_items: tuple[tuple[str, MetadataValue], ...]
    latent_items: tuple[tuple[str, MetadataValue], ...] = ()
    observational_equivalence_group: str | None = None

    def __post_init__(self) -> None:
        if not self.provider:
            raise ValueError("truth provider name must be non-empty")
        _validate_metadata_items(self.parameter_items, "parameter")
        _validate_metadata_items(self.latent_items, "latent")
        if (
            self.observational_equivalence_group is not None
            and not self.observational_equivalence_group
        ):
            raise ValueError("observational equivalence group must be non-empty")

    @property
    def parameters(self) -> dict[str, MetadataValue]:
        """Return a copy of the declared generator parameters."""

        return dict(self.parameter_items)

    @property
    def latent_metadata(self) -> dict[str, MetadataValue]:
        """Return a copy of labels intentionally absent from the density map."""

        return dict(self.latent_items)


@dataclass(frozen=True)
class AnalyticTruthSequence:
    """Time-ordered analytic truth maps for sequence-level stress tests.

    The sequence stores a declared schedule and the resulting immutable density
    maps.  It does not attach an equation of motion or imply that the schedule
    is dynamically accessible to a real condensate.
    """

    frames: tuple[AnalyticTruthInput, ...]
    times_s: tuple[float, ...]
    provider: str
    parameter_items: tuple[tuple[str, MetadataValue], ...]
    description: str

    def __post_init__(self) -> None:
        frames = tuple(self.frames)
        times_s = tuple(float(value) for value in self.times_s)
        parameter_items = tuple(self.parameter_items)
        object.__setattr__(self, "frames", frames)
        object.__setattr__(self, "times_s", times_s)
        object.__setattr__(self, "parameter_items", parameter_items)
        if not self.provider:
            raise ValueError("truth-sequence provider name must be non-empty")
        if not self.description:
            raise ValueError("truth-sequence description must be non-empty")
        if not frames or len(frames) != len(times_s):
            raise ValueError(
                "truth-sequence frames and times must be non-empty and aligned"
            )
        times = np.asarray(times_s, dtype=float)
        if np.any(~np.isfinite(times)) or np.any(np.diff(times) <= 0.0):
            raise ValueError(
                "truth-sequence times must be finite and strictly increasing"
            )
        shapes = {frame.morphology.column_density_m2.shape for frame in frames}
        names = [frame.morphology.name for frame in frames]
        if len(shapes) != 1:
            raise ValueError("truth-sequence density maps must have one common shape")
        if len(names) != len(set(names)):
            raise ValueError("truth-sequence frame names must be unique")
        if any(frame.morphology.column_density_m2.flags.writeable for frame in frames):
            raise ValueError("truth-sequence density maps must be immutable")
        _validate_metadata_items(parameter_items, "sequence parameter")

    @property
    def parameters(self) -> dict[str, MetadataValue]:
        """Return a copy of the declared sequence schedule and provenance."""

        return dict(self.parameter_items)


def _validate_metadata_items(
    items: tuple[tuple[str, MetadataValue], ...],
    kind: str,
) -> None:
    keys = [key for key, _ in items]
    if len(keys) != len(set(keys)) or any(not key for key in keys):
        raise ValueError(f"{kind} metadata keys must be non-empty and unique")
    for key, value in items:
        if isinstance(value, bool) or isinstance(value, (str, int)):
            continue
        values = (value,) if isinstance(value, float) else value
        if any(not np.isfinite(float(item)) for item in values):
            raise ValueError(f"{kind} metadata value for {key!r} must be finite")


def _metadata(**values: MetadataValue) -> tuple[tuple[str, MetadataValue], ...]:
    return tuple(values.items())


def _validated_grids(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    y = np.asarray(y_grid_m, dtype=float)
    z = np.asarray(z_grid_m, dtype=float)
    if y.ndim != 2 or y.shape != z.shape:
        raise ValueError("truth-input coordinate grids must be same-shape 2D arrays")
    if np.any(~np.isfinite(y)) or np.any(~np.isfinite(z)):
        raise ValueError("truth-input coordinate grids must be finite")
    return y * 1e6, z * 1e6


def _finite(name: str, value: float) -> float:
    checked = float(value)
    if not np.isfinite(checked):
        raise ValueError(f"{name} must be finite")
    return checked


def _positive(name: str, value: float) -> float:
    checked = _finite(name, value)
    if checked <= 0.0:
        raise ValueError(f"{name} must be positive")
    return checked


def _contrast(name: str, value: float) -> float:
    checked = _finite(name, value)
    if not 0.0 <= checked <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return checked


def _rotated_coordinates(
    y_um: NDArray[np.floating],
    z_um: NDArray[np.floating],
    *,
    centre_y_um: float,
    centre_z_um: float,
    angle_deg: float,
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    angle = np.deg2rad(_finite("angle_deg", angle_deg))
    dy = y_um - _finite("centre_y_um", centre_y_um)
    dz = z_um - _finite("centre_z_um", centre_z_um)
    longitudinal = np.cos(angle) * dy + np.sin(angle) * dz
    transverse = -np.sin(angle) * dy + np.cos(angle) * dz
    return longitudinal, transverse


def _normalised_density(
    values: ArrayLike,
    peak_column_density_m2: float,
) -> NDArray[np.floating]:
    peak = _positive("peak_column_density_m2", peak_column_density_m2)
    density = np.asarray(values, dtype=float)
    if density.ndim != 2 or np.any(~np.isfinite(density)):
        raise ValueError("analytic truth density must be a finite 2D map")
    density = np.clip(density, 0.0, None)
    maximum = float(np.max(density))
    if maximum <= np.finfo(float).tiny:
        raise ValueError("analytic truth density has no support on the supplied grid")
    result = np.array(peak * density / maximum, copy=True)
    result.setflags(write=False)
    return result


def _immutable_density(values: ArrayLike) -> NDArray[np.floating]:
    density = np.asarray(values, dtype=float)
    if density.ndim != 2 or np.any(~np.isfinite(density)):
        raise ValueError("analytic truth density must be a finite 2D map")
    if np.any(density < 0.0) or not np.any(density > 0.0):
        raise ValueError("analytic truth density must be non-negative and non-empty")
    result = np.array(density, copy=True)
    result.setflags(write=False)
    return result


def _schedule(
    name: str,
    values: tuple[float, ...] | None,
    *,
    frame_count: int,
    default: float,
    positive: bool = False,
) -> tuple[float, ...]:
    schedule = (default,) * frame_count if values is None else tuple(values)
    if len(schedule) != frame_count:
        raise ValueError(f"{name} must contain one value per sequence frame")
    checked = tuple(_finite(name, value) for value in schedule)
    if positive and any(value <= 0.0 for value in checked):
        raise ValueError(f"{name} values must be positive")
    return checked


def _rectilinear_axes_um(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
) -> tuple[
    NDArray[np.floating],
    NDArray[np.floating],
    NDArray[np.floating],
    NDArray[np.floating],
]:
    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    y_axis = np.asarray(y_um[0, :], dtype=float)
    z_axis = np.asarray(z_um[:, 0], dtype=float)
    if not np.allclose(y_um, y_axis[None, :]) or not np.allclose(
        z_um, z_axis[:, None]
    ):
        raise ValueError("sequence transforms require a rectilinear meshgrid")
    if np.any(np.diff(y_axis) <= 0.0) or np.any(np.diff(z_axis) <= 0.0):
        raise ValueError("sequence coordinate axes must be strictly increasing")
    return y_um, z_um, y_axis, z_axis


def _bilinear_sample(
    density: NDArray[np.floating],
    y_axis_um: NDArray[np.floating],
    z_axis_um: NDArray[np.floating],
    query_y_um: NDArray[np.floating],
    query_z_um: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Sample a rectilinear density map, assigning zero outside its field."""

    y_index = np.searchsorted(y_axis_um, query_y_um, side="right") - 1
    z_index = np.searchsorted(z_axis_um, query_z_um, side="right") - 1
    valid = (
        (query_y_um >= y_axis_um[0])
        & (query_y_um <= y_axis_um[-1])
        & (query_z_um >= z_axis_um[0])
        & (query_z_um <= z_axis_um[-1])
    )
    y_index = np.clip(y_index, 0, y_axis_um.size - 2)
    z_index = np.clip(z_index, 0, z_axis_um.size - 2)
    y0 = y_axis_um[y_index]
    y1 = y_axis_um[y_index + 1]
    z0 = z_axis_um[z_index]
    z1 = z_axis_um[z_index + 1]
    weight_y = (query_y_um - y0) / (y1 - y0)
    weight_z = (query_z_um - z0) / (z1 - z0)
    sampled = (
        (1.0 - weight_y) * (1.0 - weight_z) * density[z_index, y_index]
        + weight_y * (1.0 - weight_z) * density[z_index, y_index + 1]
        + (1.0 - weight_y) * weight_z * density[z_index + 1, y_index]
        + weight_y * weight_z * density[z_index + 1, y_index + 1]
    )
    return np.where(valid, sampled, 0.0)


def _truth(
    *,
    name: str,
    feature_class: str,
    description: str,
    provider: str,
    density: NDArray[np.floating],
    parameters: tuple[tuple[str, MetadataValue], ...],
) -> AnalyticTruthInput:
    return AnalyticTruthInput(
        morphology=SyntheticMorphology(name, density, description, feature_class),
        provider=provider,
        parameter_items=parameters,
    )


def rotated_asymmetric_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    longitudinal_sigma_um: float = 18.0,
    transverse_sigma_um: float = 5.0,
    angle_deg: float = 25.0,
    asymmetry: float = 0.35,
    name: str = "rotated_asymmetric_cloud",
) -> AnalyticTruthInput:
    """Return a rotated Gaussian cloud with independently controlled skew."""

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    sigma_long = _positive("longitudinal_sigma_um", longitudinal_sigma_um)
    sigma_transverse = _positive("transverse_sigma_um", transverse_sigma_um)
    skew = _finite("asymmetry", asymmetry)
    if abs(skew) > 0.95:
        raise ValueError("asymmetry magnitude must not exceed 0.95")
    longitudinal, transverse = _rotated_coordinates(
        y_um,
        z_um,
        centre_y_um=centre_y_um,
        centre_z_um=centre_z_um,
        angle_deg=angle_deg,
    )
    envelope = np.exp(
        -0.5
        * ((longitudinal / sigma_long) ** 2 + (transverse / sigma_transverse) ** 2)
    )
    values = envelope * (1.0 + skew * np.tanh(longitudinal / sigma_long))
    density = _normalised_density(values, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="rotated_asymmetry",
        description=(
            "Rotated smooth cloud with a declared longitudinal skew; analytic "
            "stress test, not an equilibrium-state prediction."
        ),
        provider="rotated_asymmetric_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            longitudinal_sigma_um=sigma_long,
            transverse_sigma_um=sigma_transverse,
            angle_deg=float(angle_deg),
            asymmetry=skew,
        ),
    )


def core_halo_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    core_sigma_y_um: float = 10.0,
    core_sigma_z_um: float = 3.0,
    halo_sigma_y_um: float = 24.0,
    halo_sigma_z_um: float = 9.0,
    halo_integrated_fraction: float = 0.20,
    core_centre_y_um: float = 0.0,
    core_centre_z_um: float = 0.0,
    halo_offset_y_um: float = 0.0,
    halo_offset_z_um: float = 0.0,
    name: str = "core_halo_cloud",
) -> AnalyticTruthInput:
    """Return a non-negative two-scale density with a declared halo share.

    Each Gaussian component is normalised on the supplied grid before mixing,
    so ``halo_integrated_fraction`` is the halo contribution to the discrete
    spatial integral before the common peak normalisation.  It is a geometric
    stress-test parameter and carries no thermodynamic interpretation.
    """

    y_um, z_um, y_axis, z_axis = _rectilinear_axes_um(y_grid_m, z_grid_m)
    core_sigma_y = _positive("core_sigma_y_um", core_sigma_y_um)
    core_sigma_z = _positive("core_sigma_z_um", core_sigma_z_um)
    halo_sigma_y = _positive("halo_sigma_y_um", halo_sigma_y_um)
    halo_sigma_z = _positive("halo_sigma_z_um", halo_sigma_z_um)
    halo_fraction = _contrast(
        "halo_integrated_fraction", halo_integrated_fraction
    )
    core_y = _finite("core_centre_y_um", core_centre_y_um)
    core_z = _finite("core_centre_z_um", core_centre_z_um)
    offset_y = _finite("halo_offset_y_um", halo_offset_y_um)
    offset_z = _finite("halo_offset_z_um", halo_offset_z_um)
    core = np.exp(
        -0.5
        * (
            ((y_um - core_y) / core_sigma_y) ** 2
            + ((z_um - core_z) / core_sigma_z) ** 2
        )
    )
    halo = np.exp(
        -0.5
        * (
            ((y_um - core_y - offset_y) / halo_sigma_y) ** 2
            + ((z_um - core_z - offset_z) / halo_sigma_z) ** 2
        )
    )

    def integrated(values: NDArray[np.floating]) -> float:
        return float(
            np.trapezoid(np.trapezoid(values, y_axis, axis=1), z_axis, axis=0)
        )

    values = (1.0 - halo_fraction) * core / integrated(core)
    values += halo_fraction * halo / integrated(halo)
    density = _normalised_density(values, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="two_scale_core_halo",
        description=(
            "Controlled two-scale core and halo density with an explicit "
            "integrated halo share and displacement; no population or "
            "thermodynamic interpretation is attached."
        ),
        provider="core_halo_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            core_sigma_y_um=core_sigma_y,
            core_sigma_z_um=core_sigma_z,
            halo_sigma_y_um=halo_sigma_y,
            halo_sigma_z_um=halo_sigma_z,
            halo_integrated_fraction=halo_fraction,
            core_centre_y_um=core_y,
            core_centre_z_um=core_z,
            halo_offset_y_um=offset_y,
            halo_offset_z_um=offset_z,
        ),
    )


def transverse_modulated_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    longitudinal_sigma_um: float = 20.0,
    transverse_sigma_um: float = 7.0,
    modulation_period_um: float = 5.0,
    modulation_contrast: float = 0.25,
    modulation_phase_rad: float = 0.0,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    name: str = "transverse_modulated_cloud",
) -> AnalyticTruthInput:
    """Return a smooth cloud carrying modulation along the transverse axis."""

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    sigma_long = _positive("longitudinal_sigma_um", longitudinal_sigma_um)
    sigma_transverse = _positive("transverse_sigma_um", transverse_sigma_um)
    period = _positive("modulation_period_um", modulation_period_um)
    contrast = _contrast("modulation_contrast", modulation_contrast)
    phase = _finite("modulation_phase_rad", modulation_phase_rad)
    dy = y_um - _finite("centre_y_um", centre_y_um)
    dz = z_um - _finite("centre_z_um", centre_z_um)
    envelope = np.exp(
        -0.5 * ((dy / sigma_long) ** 2 + (dz / sigma_transverse) ** 2)
    )
    modulation = 1.0 + contrast * np.cos(2.0 * np.pi * dz / period + phase)
    density = _normalised_density(envelope * modulation, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="transverse_modulation",
        description=(
            "Smooth envelope with controlled transverse period, phase and "
            "contrast; not a roton or supersolid-state claim."
        ),
        provider="transverse_modulated_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            longitudinal_sigma_um=sigma_long,
            transverse_sigma_um=sigma_transverse,
            modulation_period_um=period,
            modulation_contrast=contrast,
            modulation_phase_rad=phase,
        ),
    )


def gaussian_peak_chain(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    peak_count: int = 3,
    separation_um: float = 9.0,
    longitudinal_peak_sigma_um: float = 2.2,
    transverse_peak_sigma_um: float = 2.0,
    peak_contrast: float = 0.85,
    angle_deg: float = 0.0,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    relative_peak_weights: tuple[float, ...] | None = None,
    name: str | None = None,
) -> AnalyticTruthInput:
    """Return a one- to five-peak chain with independent scale and contrast."""

    if isinstance(peak_count, bool) or int(peak_count) != peak_count:
        raise ValueError("peak_count must be an integer")
    count = int(peak_count)
    if not 1 <= count <= 5:
        raise ValueError("peak_count must lie between one and five")
    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    separation = _positive("separation_um", separation_um)
    sigma_long = _positive(
        "longitudinal_peak_sigma_um", longitudinal_peak_sigma_um
    )
    sigma_transverse = _positive(
        "transverse_peak_sigma_um", transverse_peak_sigma_um
    )
    contrast = _contrast("peak_contrast", peak_contrast)
    if relative_peak_weights is None:
        weights = np.ones(count, dtype=float)
    else:
        weights = np.asarray(relative_peak_weights, dtype=float)
        if weights.shape != (count,):
            raise ValueError("relative_peak_weights must match peak_count")
        if np.any(~np.isfinite(weights)) or np.any(weights < 0.0) or not np.any(weights > 0.0):
            raise ValueError("relative peak weights must be finite and non-negative")
    longitudinal, transverse = _rotated_coordinates(
        y_um,
        z_um,
        centre_y_um=centre_y_um,
        centre_z_um=centre_z_um,
        angle_deg=angle_deg,
    )
    positions = separation * (np.arange(count, dtype=float) - 0.5 * (count - 1))
    peaks = np.zeros_like(longitudinal)
    for position, weight in zip(positions, weights, strict=True):
        peaks += weight * np.exp(
            -0.5
            * (
                ((longitudinal - position) / sigma_long) ** 2
                + (transverse / sigma_transverse) ** 2
            )
        )
    chain_sigma = max(sigma_long, 0.5 * separation * max(count - 1, 1) + sigma_long)
    pedestal = np.exp(
        -0.5
        * ((longitudinal / chain_sigma) ** 2 + (transverse / sigma_transverse) ** 2)
    )
    density = _normalised_density(
        contrast * peaks + (1.0 - contrast) * pedestal,
        peak_column_density_m2,
    )
    label = name or f"gaussian_{count}_peak_chain"
    return _truth(
        name=label,
        feature_class="peak_chain",
        description=(
            f"Controlled analytic {count}-peak chain with independent peak "
            "width, separation and pedestal contrast."
        ),
        provider="gaussian_peak_chain",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            peak_count=count,
            separation_um=separation,
            longitudinal_peak_sigma_um=sigma_long,
            transverse_peak_sigma_um=sigma_transverse,
            peak_contrast=contrast,
            angle_deg=float(angle_deg),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            relative_peak_weights=tuple(float(value) for value in weights),
        ),
    )


def central_hole_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    longitudinal_sigma_um: float = 15.0,
    transverse_sigma_um: float = 6.0,
    hole_sigma_um: float = 2.5,
    hole_contrast: float = 1.0,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    angle_deg: float = 0.0,
    name: str = "central_hole_cloud",
) -> AnalyticTruthInput:
    """Return an elliptical cloud with a controlled central hole or ring."""

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    sigma_long = _positive("longitudinal_sigma_um", longitudinal_sigma_um)
    sigma_transverse = _positive("transverse_sigma_um", transverse_sigma_um)
    hole_sigma = _positive("hole_sigma_um", hole_sigma_um)
    contrast = _contrast("hole_contrast", hole_contrast)
    longitudinal, transverse = _rotated_coordinates(
        y_um,
        z_um,
        centre_y_um=centre_y_um,
        centre_z_um=centre_z_um,
        angle_deg=angle_deg,
    )
    envelope = np.exp(
        -0.5
        * ((longitudinal / sigma_long) ** 2 + (transverse / sigma_transverse) ** 2)
    )
    radial_distance = np.sqrt(longitudinal**2 + transverse**2)
    hole = 1.0 - contrast * np.exp(-0.5 * (radial_distance / hole_sigma) ** 2)
    density = _normalised_density(envelope * hole, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="ring_or_hole",
        description=(
            "Elliptical density envelope with a controlled central depletion; "
            "the map does not encode circulation or a condensate phase winding."
        ),
        provider="central_hole_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            longitudinal_sigma_um=sigma_long,
            transverse_sigma_um=sigma_transverse,
            hole_sigma_um=hole_sigma,
            hole_contrast=contrast,
            angle_deg=float(angle_deg),
        ),
    )


def curved_notched_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    longitudinal_sigma_um: float = 20.0,
    transverse_sigma_um: float = 2.5,
    curvature_displacement_um: float = 3.0,
    notch_centre_y_um: float = 3.0,
    notch_sigma_um: float = 2.5,
    notch_contrast: float = 0.7,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    name: str = "curved_notched_cloud",
) -> AnalyticTruthInput:
    """Return a curved density ridge with a local longitudinal notch."""

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    sigma_long = _positive("longitudinal_sigma_um", longitudinal_sigma_um)
    sigma_transverse = _positive("transverse_sigma_um", transverse_sigma_um)
    curvature = _finite("curvature_displacement_um", curvature_displacement_um)
    notch_centre = _finite("notch_centre_y_um", notch_centre_y_um)
    notch_sigma = _positive("notch_sigma_um", notch_sigma_um)
    contrast = _contrast("notch_contrast", notch_contrast)
    dy = y_um - _finite("centre_y_um", centre_y_um)
    centreline = _finite("centre_z_um", centre_z_um) + curvature * (
        dy / sigma_long
    ) ** 2
    transverse_offset = z_um - centreline
    ridge = np.exp(
        -0.5 * ((dy / sigma_long) ** 2 + (transverse_offset / sigma_transverse) ** 2)
    )
    notch = 1.0 - contrast * np.exp(
        -0.5 * ((y_um - notch_centre) / notch_sigma) ** 2
    )
    density = _normalised_density(ridge * notch, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="curved_local_defect",
        description=(
            "Curved analytic density ridge with a controlled local notch; no "
            "soliton phase or dynamical interpretation is attached."
        ),
        provider="curved_notched_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            longitudinal_sigma_um=sigma_long,
            transverse_sigma_um=sigma_transverse,
            curvature_displacement_um=curvature,
            notch_centre_y_um=notch_centre,
            notch_sigma_um=notch_sigma,
            notch_contrast=contrast,
        ),
    )


def expanded_shifted_cloud(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    centre_y_um: float = 24.0,
    centre_z_um: float = 5.0,
    longitudinal_sigma_um: float = 24.0,
    transverse_sigma_um: float = 9.0,
    angle_deg: float = 12.0,
    super_gaussian_order: float = 4.0,
    name: str = "expanded_shifted_support_challenge",
) -> AnalyticTruthInput:
    """Return a broad shifted cloud intended to challenge a fixed fit support."""

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    sigma_long = _positive("longitudinal_sigma_um", longitudinal_sigma_um)
    sigma_transverse = _positive("transverse_sigma_um", transverse_sigma_um)
    order = _positive("super_gaussian_order", super_gaussian_order)
    longitudinal, transverse = _rotated_coordinates(
        y_um,
        z_um,
        centre_y_um=centre_y_um,
        centre_z_um=centre_z_um,
        angle_deg=angle_deg,
    )
    values = np.exp(
        -0.5
        * (
            np.abs(longitudinal / sigma_long) ** order
            + np.abs(transverse / sigma_transverse) ** order
        )
    )
    density = _normalised_density(values, peak_column_density_m2)
    return _truth(
        name=name,
        feature_class="support_challenge",
        description=(
            "Broad shifted super-Gaussian cloud intended to expose fixed-support "
            "and boundary-prior bias."
        ),
        provider="expanded_shifted_cloud",
        density=density,
        parameters=_metadata(
            peak_column_density_m2=float(peak_column_density_m2),
            centre_y_um=float(centre_y_um),
            centre_z_um=float(centre_z_um),
            longitudinal_sigma_um=sigma_long,
            transverse_sigma_um=sigma_transverse,
            angle_deg=float(angle_deg),
            super_gaussian_order=order,
        ),
    )


def controlled_unknown_sequence(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    base: AnalyticTruthInput | ArrayLike,
    times_s: tuple[float, ...],
    amplitude_scales: tuple[float, ...] | None = None,
    translation_y_um: tuple[float, ...] | None = None,
    translation_z_um: tuple[float, ...] | None = None,
    scale_y: tuple[float, ...] | None = None,
    scale_z: tuple[float, ...] | None = None,
    transform_centre_y_um: float = 0.0,
    transform_centre_z_um: float = 0.0,
    sudden_peak_frame_index: int | None = None,
    sudden_peak_relative_amplitude: float = 0.25,
    sudden_peak_centre_y_um: float = 8.0,
    sudden_peak_centre_z_um: float = 0.0,
    sudden_peak_sigma_y_um: float = 2.0,
    sudden_peak_sigma_z_um: float = 2.0,
    name: str = "controlled_unknown_sequence",
) -> AnalyticTruthSequence:
    """Build a controlled time sequence without assuming condensate dynamics.

    The declared schedules independently change amplitude, translation and
    spatial scale.  An optional positive Gaussian feature appears at the named
    frame and persists thereafter.  These operations define synthetic truth
    for reconstruction tests only; no forward model, fitter or equation of
    motion is called.
    """

    y_um, z_um, y_axis, z_axis = _rectilinear_axes_um(y_grid_m, z_grid_m)
    times = tuple(_finite("times_s", value) for value in times_s)
    if not times or any(later <= earlier for earlier, later in zip(times, times[1:])):
        raise ValueError("times_s must be non-empty and strictly increasing")
    frame_count = len(times)
    amplitudes = _schedule(
        "amplitude_scales",
        amplitude_scales,
        frame_count=frame_count,
        default=1.0,
        positive=True,
    )
    shifts_y = _schedule(
        "translation_y_um",
        translation_y_um,
        frame_count=frame_count,
        default=0.0,
    )
    shifts_z = _schedule(
        "translation_z_um",
        translation_z_um,
        frame_count=frame_count,
        default=0.0,
    )
    scales_y = _schedule(
        "scale_y", scale_y, frame_count=frame_count, default=1.0, positive=True
    )
    scales_z = _schedule(
        "scale_z", scale_z, frame_count=frame_count, default=1.0, positive=True
    )
    centre_y = _finite("transform_centre_y_um", transform_centre_y_um)
    centre_z = _finite("transform_centre_z_um", transform_centre_z_um)

    if isinstance(base, AnalyticTruthInput):
        source = base.morphology.column_density_m2
        source_provider = base.provider
        source_name = base.morphology.name
    else:
        source = _immutable_density(base)
        source_provider = "array_input"
        source_name = "unnamed_density_map"
    if source.shape != y_um.shape:
        raise ValueError("base density must have the same shape as the sequence grid")

    if sudden_peak_frame_index is not None:
        if (
            isinstance(sudden_peak_frame_index, bool)
            or int(sudden_peak_frame_index) != sudden_peak_frame_index
            or not 0 <= int(sudden_peak_frame_index) < frame_count
        ):
            raise ValueError("sudden_peak_frame_index must select a sequence frame")
        sudden_index = int(sudden_peak_frame_index)
        sudden_amplitude = _positive(
            "sudden_peak_relative_amplitude", sudden_peak_relative_amplitude
        )
        sudden_sigma_y = _positive("sudden_peak_sigma_y_um", sudden_peak_sigma_y_um)
        sudden_sigma_z = _positive("sudden_peak_sigma_z_um", sudden_peak_sigma_z_um)
    else:
        sudden_index = None
        sudden_amplitude = _finite(
            "sudden_peak_relative_amplitude", sudden_peak_relative_amplitude
        )
        sudden_sigma_y = _positive("sudden_peak_sigma_y_um", sudden_peak_sigma_y_um)
        sudden_sigma_z = _positive("sudden_peak_sigma_z_um", sudden_peak_sigma_z_um)
    sudden_y = _finite("sudden_peak_centre_y_um", sudden_peak_centre_y_um)
    sudden_z = _finite("sudden_peak_centre_z_um", sudden_peak_centre_z_um)
    peak_template = np.exp(
        -0.5
        * (
            ((y_um - sudden_y) / sudden_sigma_y) ** 2
            + ((z_um - sudden_z) / sudden_sigma_z) ** 2
        )
    )
    source_peak = float(np.max(source))

    frames: list[AnalyticTruthInput] = []
    for index, (time_s, amplitude, shift_y, shift_z, stretch_y, stretch_z) in enumerate(
        zip(
            times,
            amplitudes,
            shifts_y,
            shifts_z,
            scales_y,
            scales_z,
            strict=True,
        )
    ):
        query_y = centre_y + (y_um - centre_y - shift_y) / stretch_y
        query_z = centre_z + (z_um - centre_z - shift_z) / stretch_z
        values = amplitude * _bilinear_sample(
            source, y_axis, z_axis, query_y, query_z
        )
        peak_active = sudden_index is not None and index >= sudden_index
        if peak_active:
            values = values + sudden_amplitude * source_peak * peak_template
        density = _immutable_density(values)
        frames.append(
            _truth(
                name=f"{name}_frame_{index:03d}",
                feature_class="controlled_unknown_sequence_frame",
                description=(
                    "Controlled analytic sequence frame generated by declared "
                    "amplitude and coordinate schedules, without a dynamical law."
                ),
                provider="controlled_unknown_sequence",
                density=density,
                parameters=_metadata(
                    source_provider=source_provider,
                    source_name=source_name,
                    frame_index=index,
                    time_s=time_s,
                    amplitude_scale=amplitude,
                    translation_y_um=shift_y,
                    translation_z_um=shift_z,
                    scale_y=stretch_y,
                    scale_z=stretch_z,
                    sudden_peak_active=peak_active,
                ),
            )
        )

    return AnalyticTruthSequence(
        frames=tuple(frames),
        times_s=times,
        provider="controlled_unknown_sequence",
        description=(
            "Controlled amplitude, translation and scale schedules with an "
            "optional abrupt feature; no condensate dynamics are assumed."
        ),
        parameter_items=_metadata(
            source_provider=source_provider,
            source_name=source_name,
            amplitude_scales=amplitudes,
            translation_y_um=shifts_y,
            translation_z_um=shifts_z,
            scale_y=scales_y,
            scale_z=scales_z,
            transform_centre_y_um=centre_y,
            transform_centre_z_um=centre_z,
            sudden_peak_frame_index=(
                sudden_index if sudden_index is not None else "none"
            ),
            sudden_peak_relative_amplitude=sudden_amplitude,
            sudden_peak_centre_y_um=sudden_y,
            sudden_peak_centre_z_um=sudden_z,
            sudden_peak_sigma_y_um=sudden_sigma_y,
            sudden_peak_sigma_z_um=sudden_sigma_z,
            dynamics_model="none",
        ),
    )


def latent_phase_impossibility_pair(
    base: AnalyticTruthInput,
    *,
    phase_labels: tuple[str, str] = ("phase_label_A", "phase_label_B"),
    equivalence_group: str | None = None,
) -> tuple[AnalyticTruthInput, AnalyticTruthInput]:
    """Attach two different latent phase labels to one identical density map.

    The returned inputs deliberately share the same immutable column-density
    array.  A density-only forward model must therefore produce identical
    observations for both.  The labels are metadata hooks for an explicit
    non-identifiability test, not simulated condensate phases.
    """

    if len(phase_labels) != 2 or any(not label for label in phase_labels):
        raise ValueError("phase_labels must contain two non-empty labels")
    if phase_labels[0] == phase_labels[1]:
        raise ValueError("impossibility-pair phase labels must differ")
    group = equivalence_group or f"{base.morphology.name}__same_density_phase_pair"
    density = base.morphology.column_density_m2
    pair: list[AnalyticTruthInput] = []
    for index, label in enumerate(phase_labels):
        morphology = SyntheticMorphology(
            name=f"{base.morphology.name}__latent_{index}",
            column_density_m2=density,
            description=(
                f"{base.morphology.description} Identical density with latent "
                f"phase label {label!r}."
            ),
            feature_class=base.morphology.feature_class,
        )
        pair.append(
            AnalyticTruthInput(
                morphology=morphology,
                provider=base.provider,
                parameter_items=base.parameter_items,
                latent_items=_metadata(
                    latent_phase_label=label,
                    density_observability="identical_column_density",
                ),
                observational_equivalence_group=group,
            )
        )
    return pair[0], pair[1]


def build_analytic_truth_catalogue(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
) -> tuple[AnalyticTruthInput, ...]:
    """Return one representative input from each analytic stress-test family."""

    return (
        rotated_asymmetric_cloud(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
        ),
        transverse_modulated_cloud(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
        ),
        gaussian_peak_chain(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
            peak_count=3,
        ),
        central_hole_cloud(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
        ),
        curved_notched_cloud(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
        ),
        expanded_shifted_cloud(
            y_grid_m,
            z_grid_m,
            peak_column_density_m2=peak_column_density_m2,
        ),
    )
