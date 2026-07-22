"""Grid construction and convergence checks for inverse studies.

The canonical forward calculation can use either the historical integer block
average or an exact area-weighted physical-camera sampler.  Repeating a
nonlinear inverse fit on the full object grid is unnecessarily expensive: the
optical aperture and camera sampling remove most of the additional numerical
degrees of freedom.  This module constructs smaller grids whose camera pixels
coincide with the canonical physical coordinates and quantifies any remaining
readout difference.  A reduced grid is therefore an explicitly checked
numerical approximation, not a change to the optical model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .contracts import ReconstructionGrid


@dataclass(frozen=True)
class CameraAlignedReduction:
    """Provenance for a reduced grid aligned to canonical camera pixels."""

    canonical_ngrid: int
    canonical_field_of_view_m: float
    canonical_bin_size: int
    canonical_trimmed_cells: int
    camera_pixel_count: int
    reduced_ngrid: int
    reduced_field_of_view_m: float
    reduced_bin_size: int
    reduced_coordinate_shift_m: float
    maximum_camera_coordinate_mismatch_m: float


@dataclass(frozen=True)
class PhysicalCameraAlignedReduction:
    """Provenance for a reduced grid using the canonical physical camera."""

    canonical_ngrid: int
    canonical_field_of_view_m: float
    camera_pixel_size_m: float
    camera_output_shape: tuple[int, int]
    reduced_ngrid: int
    reduced_field_of_view_m: float
    maximum_camera_coordinate_mismatch_m: float


@dataclass(frozen=True)
class FaradayGridAgreement:
    """Signal-level agreement between canonical and reduced-grid readouts."""

    dual_port_signal_relative_l2_error: float
    dual_port_atom_dependent_channels_relative_l2_error: float
    dark_field_relative_l2_error: float
    dual_port_peak_relative_error: float
    dark_field_peak_relative_error: float


def _uniform_object_arrays(
    *,
    ngrid: int,
    field_of_view_m: float,
    numerical_aperture: float,
    wavelength_m: float,
    coordinate_shift_m: float = 0.0,
) -> tuple[
    NDArray[np.floating],
    NDArray[np.floating],
    NDArray[np.floating],
]:
    """Build the shared object coordinates and circular pupil."""

    if ngrid <= 0:
        raise ValueError("ngrid must be positive")
    if not np.isfinite(field_of_view_m) or field_of_view_m <= 0:
        raise ValueError("field of view must be finite and positive")
    if not np.isfinite(numerical_aperture) or numerical_aperture <= 0:
        raise ValueError("numerical aperture must be finite and positive")
    if not np.isfinite(wavelength_m) or wavelength_m <= 0:
        raise ValueError("wavelength must be finite and positive")
    if not np.isfinite(coordinate_shift_m):
        raise ValueError("coordinate shift must be finite")

    spacing_m = field_of_view_m / ngrid
    coordinate_axis_m = (
        (np.arange(ngrid, dtype=float) - ngrid // 2) * spacing_m
        + coordinate_shift_m
    )
    y_grid_m, z_grid_m = np.meshgrid(coordinate_axis_m, coordinate_axis_m)
    frequency_axis_m_inv = np.fft.fftfreq(ngrid, d=spacing_m)
    frequency_y_m_inv, frequency_z_m_inv = np.meshgrid(
        frequency_axis_m_inv,
        frequency_axis_m_inv,
    )
    pupil = (
        np.sqrt(frequency_y_m_inv**2 + frequency_z_m_inv**2)
        <= numerical_aperture / wavelength_m
    ).astype(float)
    return y_grid_m, z_grid_m, pupil


def build_uniform_reconstruction_grid(
    *,
    ngrid: int,
    field_of_view_m: float,
    bin_size: int,
    numerical_aperture: float,
    wavelength_m: float,
    coordinate_shift_m: float = 0.0,
    roi_half_width_y_um: float | None = None,
    roi_half_width_z_um: float | None = None,
) -> ReconstructionGrid:
    """Build the periodic object grid, circular pupil and fixed camera ROI."""

    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    if (roi_half_width_y_um is None) != (roi_half_width_z_um is None):
        raise ValueError("both ROI half-widths must be supplied together")
    if roi_half_width_y_um is not None and (
        not np.isfinite(roi_half_width_y_um)
        or not np.isfinite(roi_half_width_z_um)
        or roi_half_width_y_um <= 0
        or roi_half_width_z_um <= 0
    ):
        raise ValueError("ROI half-widths must be finite and positive")

    y_grid_m, z_grid_m, pupil = _uniform_object_arrays(
        ngrid=ngrid,
        field_of_view_m=field_of_view_m,
        numerical_aperture=numerical_aperture,
        wavelength_m=wavelength_m,
        coordinate_shift_m=coordinate_shift_m,
    )
    coordinate_axis_m = y_grid_m[0, :]

    camera_rows = ngrid // bin_size
    camera_shape = (camera_rows, camera_rows)
    if roi_half_width_y_um is None:
        roi_mask = np.ones(camera_shape, dtype=bool)
    else:
        trimmed = camera_rows * bin_size
        camera_axis_um = coordinate_axis_m[:trimmed].reshape(camera_rows, bin_size).mean(axis=1) * 1e6
        camera_y_um, camera_z_um = np.meshgrid(camera_axis_um, camera_axis_um)
        roi_mask = (
            (np.abs(camera_y_um) <= float(roi_half_width_y_um))
            & (np.abs(camera_z_um) <= float(roi_half_width_z_um))
        )

    return ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=bin_size,
        roi_mask=roi_mask,
    )


def build_uniform_physical_camera_grid(
    *,
    ngrid: int,
    field_of_view_m: float,
    camera_pixel_size_m: float,
    camera_output_shape: tuple[int, int],
    numerical_aperture: float,
    wavelength_m: float,
    roi_half_width_y_um: float | None = None,
    roi_half_width_z_um: float | None = None,
) -> ReconstructionGrid:
    """Build an object grid sampled by centred physical camera pixels."""

    if len(camera_output_shape) != 2 or any(value <= 0 for value in camera_output_shape):
        raise ValueError("camera_output_shape must contain two positive dimensions")
    if not np.isfinite(camera_pixel_size_m) or camera_pixel_size_m <= 0:
        raise ValueError("camera_pixel_size_m must be finite and positive")
    if (roi_half_width_y_um is None) != (roi_half_width_z_um is None):
        raise ValueError("both ROI half-widths must be supplied together")
    if roi_half_width_y_um is not None and (
        not np.isfinite(roi_half_width_y_um)
        or not np.isfinite(roi_half_width_z_um)
        or roi_half_width_y_um <= 0
        or roi_half_width_z_um <= 0
    ):
        raise ValueError("ROI half-widths must be finite and positive")

    y_grid_m, z_grid_m, pupil = _uniform_object_arrays(
        ngrid=ngrid,
        field_of_view_m=field_of_view_m,
        numerical_aperture=numerical_aperture,
        wavelength_m=wavelength_m,
    )
    if roi_half_width_y_um is None:
        roi_mask = np.ones(camera_output_shape, dtype=bool)
    else:
        camera_y_um = (
            np.arange(camera_output_shape[1], dtype=float)
            - (camera_output_shape[1] - 1) / 2
        ) * camera_pixel_size_m * 1e6
        camera_z_um = (
            np.arange(camera_output_shape[0], dtype=float)
            - (camera_output_shape[0] - 1) / 2
        ) * camera_pixel_size_m * 1e6
        camera_y_grid_um, camera_z_grid_um = np.meshgrid(camera_y_um, camera_z_um)
        roi_mask = (
            (np.abs(camera_y_grid_um) <= float(roi_half_width_y_um))
            & (np.abs(camera_z_grid_um) <= float(roi_half_width_z_um))
        )
    return ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=None,
        roi_mask=roi_mask,
        camera_pixel_size_m=camera_pixel_size_m,
        camera_output_shape=camera_output_shape,
    )


def build_camera_aligned_reduced_grid(
    *,
    canonical_ngrid: int,
    canonical_field_of_view_m: float,
    canonical_bin_size: int,
    reduced_bin_size: int,
    numerical_aperture: float,
    wavelength_m: float,
    roi_half_width_y_um: float | None = None,
    roi_half_width_z_um: float | None = None,
) -> tuple[ReconstructionGrid, CameraAlignedReduction]:
    """Return a smaller grid with exactly coincident binned camera coordinates.

    The canonical camera discards trailing cells that do not fill a complete
    bin.  The reduced field of view represents exactly those retained cells.
    A sub-grid coordinate shift then aligns every reduced camera-pixel centre
    to the canonical centre.  The Fourier-domain sampling changes slightly and
    must still be checked at the readout level with :func:`faraday_grid_agreement`.
    """

    if canonical_ngrid <= 0 or canonical_bin_size <= 0 or reduced_bin_size <= 0:
        raise ValueError("grid sizes and bin sizes must be positive")
    canonical_trimmed_cells = (
        canonical_ngrid // canonical_bin_size
    ) * canonical_bin_size
    if canonical_trimmed_cells == 0:
        raise ValueError("canonical bin size is larger than its grid")
    camera_pixel_count = canonical_trimmed_cells // canonical_bin_size
    reduced_ngrid = camera_pixel_count * reduced_bin_size
    canonical_spacing_m = canonical_field_of_view_m / canonical_ngrid
    reduced_field_of_view_m = canonical_spacing_m * canonical_trimmed_cells
    reduced_spacing_m = reduced_field_of_view_m / reduced_ngrid

    canonical_axis_m = (
        np.arange(canonical_ngrid, dtype=float) - canonical_ngrid // 2
    ) * canonical_spacing_m
    canonical_first_camera_m = float(np.mean(canonical_axis_m[:canonical_bin_size]))
    unshifted_reduced_axis_m = (
        np.arange(reduced_ngrid, dtype=float) - reduced_ngrid // 2
    ) * reduced_spacing_m
    reduced_coordinate_shift_m = canonical_first_camera_m - float(
        np.mean(unshifted_reduced_axis_m[:reduced_bin_size])
    )

    reduced_grid = build_uniform_reconstruction_grid(
        ngrid=reduced_ngrid,
        field_of_view_m=reduced_field_of_view_m,
        bin_size=reduced_bin_size,
        numerical_aperture=numerical_aperture,
        wavelength_m=wavelength_m,
        coordinate_shift_m=reduced_coordinate_shift_m,
        roi_half_width_y_um=roi_half_width_y_um,
        roi_half_width_z_um=roi_half_width_z_um,
    )
    canonical_camera_axis_m = canonical_axis_m[:canonical_trimmed_cells].reshape(
        camera_pixel_count,
        canonical_bin_size,
    ).mean(axis=1)
    maximum_mismatch_m = float(
        np.max(np.abs(canonical_camera_axis_m - reduced_grid.camera_y_um * 1e-6))
    )
    return reduced_grid, CameraAlignedReduction(
        canonical_ngrid=int(canonical_ngrid),
        canonical_field_of_view_m=float(canonical_field_of_view_m),
        canonical_bin_size=int(canonical_bin_size),
        canonical_trimmed_cells=int(canonical_trimmed_cells),
        camera_pixel_count=int(camera_pixel_count),
        reduced_ngrid=int(reduced_ngrid),
        reduced_field_of_view_m=float(reduced_field_of_view_m),
        reduced_bin_size=int(reduced_bin_size),
        reduced_coordinate_shift_m=float(reduced_coordinate_shift_m),
        maximum_camera_coordinate_mismatch_m=maximum_mismatch_m,
    )


def build_physical_camera_aligned_reduced_grid(
    *,
    canonical_ngrid: int,
    canonical_field_of_view_m: float,
    camera_pixel_size_m: float,
    camera_output_shape: tuple[int, int],
    reduced_ngrid: int,
    numerical_aperture: float,
    wavelength_m: float,
    roi_half_width_y_um: float | None = None,
    roi_half_width_z_um: float | None = None,
) -> tuple[ReconstructionGrid, PhysicalCameraAlignedReduction]:
    """Build a smaller object grid with the same physical camera sampler.

    The reduced calculation keeps the canonical field of view and integrates
    both grids over identical physical camera pixels.  It therefore avoids
    replacing a non-integer camera pitch by an approximate block average.
    Signal-level convergence still has to be checked for every declared
    reduction.
    """

    if len(camera_output_shape) != 2 or camera_output_shape[0] != camera_output_shape[1]:
        raise ValueError("the current Fourier grid requires a square camera output")
    if reduced_ngrid <= 0:
        raise ValueError("reduced_ngrid must be positive")
    reduced_grid = build_uniform_physical_camera_grid(
        ngrid=reduced_ngrid,
        field_of_view_m=canonical_field_of_view_m,
        camera_pixel_size_m=camera_pixel_size_m,
        camera_output_shape=camera_output_shape,
        numerical_aperture=numerical_aperture,
        wavelength_m=wavelength_m,
        roi_half_width_y_um=roi_half_width_y_um,
        roi_half_width_z_um=roi_half_width_z_um,
    )
    camera_pixel_count = int(camera_output_shape[0])
    target_camera_axis_m = (
        np.arange(camera_pixel_count, dtype=float) - (camera_pixel_count - 1) / 2
    ) * camera_pixel_size_m
    maximum_mismatch_m = float(
        np.max(np.abs(target_camera_axis_m - reduced_grid.camera_y_um * 1e-6))
    )
    return reduced_grid, PhysicalCameraAlignedReduction(
        canonical_ngrid=int(canonical_ngrid),
        canonical_field_of_view_m=float(canonical_field_of_view_m),
        camera_pixel_size_m=float(camera_pixel_size_m),
        camera_output_shape=(int(camera_output_shape[0]), int(camera_output_shape[1])),
        reduced_ngrid=int(reduced_ngrid),
        reduced_field_of_view_m=float(canonical_field_of_view_m),
        maximum_camera_coordinate_mismatch_m=maximum_mismatch_m,
    )


def _relative_l2(reference: NDArray[np.floating], candidate: NDArray[np.floating]) -> float:
    denominator = max(float(np.linalg.norm(reference)), np.finfo(float).eps)
    return float(np.linalg.norm(candidate - reference) / denominator)


def faraday_grid_agreement(
    *,
    canonical_h: ArrayLike,
    canonical_v: ArrayLike,
    reduced_h: ArrayLike,
    reduced_v: ArrayLike,
    canonical_dark: ArrayLike,
    reduced_dark: ArrayLike,
) -> FaradayGridAgreement:
    """Compare the atom-dependent Faraday signals on two camera grids."""

    arrays = [
        np.asarray(value, dtype=float)
        for value in (
            canonical_h,
            canonical_v,
            reduced_h,
            reduced_v,
            canonical_dark,
            reduced_dark,
        )
    ]
    if any(array.ndim != 2 for array in arrays):
        raise ValueError("Faraday grid comparison expects two-dimensional camera maps")
    if len({array.shape for array in arrays}) != 1:
        raise ValueError("canonical and reduced camera maps must have identical shapes")
    if any(np.any(~np.isfinite(array)) for array in arrays):
        raise ValueError("Faraday grid comparison maps must be finite")
    canonical_h_array, canonical_v_array, reduced_h_array, reduced_v_array, canonical_dark_array, reduced_dark_array = arrays
    canonical_total = canonical_h_array + canonical_v_array
    reduced_total = reduced_h_array + reduced_v_array
    canonical_signal = np.divide(
        canonical_h_array - canonical_v_array,
        canonical_total,
        out=np.zeros_like(canonical_total),
        where=np.abs(canonical_total) > np.finfo(float).eps,
    )
    reduced_signal = np.divide(
        reduced_h_array - reduced_v_array,
        reduced_total,
        out=np.zeros_like(reduced_total),
        where=np.abs(reduced_total) > np.finfo(float).eps,
    )
    canonical_atom_ports = np.concatenate(
        [(canonical_h_array - 0.5).ravel(), (canonical_v_array - 0.5).ravel()]
    )
    reduced_atom_ports = np.concatenate(
        [(reduced_h_array - 0.5).ravel(), (reduced_v_array - 0.5).ravel()]
    )
    canonical_signal_peak = max(float(np.max(np.abs(canonical_signal))), np.finfo(float).eps)
    canonical_dark_peak = max(float(np.max(canonical_dark_array)), np.finfo(float).eps)
    return FaradayGridAgreement(
        dual_port_signal_relative_l2_error=_relative_l2(canonical_signal, reduced_signal),
        dual_port_atom_dependent_channels_relative_l2_error=_relative_l2(
            canonical_atom_ports,
            reduced_atom_ports,
        ),
        dark_field_relative_l2_error=_relative_l2(canonical_dark_array, reduced_dark_array),
        dual_port_peak_relative_error=float(
            np.max(np.abs(reduced_signal)) / canonical_signal_peak - 1.0
        ),
        dark_field_peak_relative_error=float(
            np.max(reduced_dark_array) / canonical_dark_peak - 1.0
        ),
    )
