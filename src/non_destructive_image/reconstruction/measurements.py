"""Raw-channel Faraday measurement operators for inverse fitting."""

from __future__ import annotations

from collections.abc import Iterable
from typing import ClassVar, Protocol, runtime_checkable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import fft as scipy_fft

from ..fourier import propagate_scattered_field
from .contracts import DetectorContract, FaradayResponseContract, ReconstructionGrid
from .noise import simulate_poisson_gaussian_counts
from .object_models import (
    DifferentiableColumnDensityModel,
    smooth_tf_column_density,
    smooth_tf_density_and_internal_jacobian,
)
from .parameters import SmoothTFParameters, from_internal


@runtime_checkable
class DifferentiableDensityMeasurement(Protocol):
    """Raw-count instrument contract consumed by the generic density fit.

    Faraday readouts implement this contract now.  Future PCI, DGI or
    absorption operators can implement the same methods without changing the
    object models or optimisation layer.
    """

    grid: ReconstructionGrid

    @property
    def read_noise_electrons(self) -> float:
        ...

    def flatten_observed(self, *channels: ArrayLike) -> NDArray[np.floating]:
        ...

    def expected_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], ...]:
        ...

    def expected_vector_and_jacobian_model(
        self,
        model: DifferentiableColumnDensityModel,
        parameter_vector: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        ...


class _FaradayMeasurementBase:
    """Common finite-aperture propagation for the two Faraday readouts."""

    channel_names: ClassVar[tuple[str, ...]]

    def __init__(
        self,
        *,
        grid: ReconstructionGrid,
        detector: DetectorContract,
        response: FaradayResponseContract,
        jacobian_batch_size: int = 4,
    ) -> None:
        if jacobian_batch_size <= 0:
            raise ValueError("Jacobian batch size must be positive")
        self.grid = grid
        self.detector = detector
        self.response = response
        self.jacobian_batch_size = int(jacobian_batch_size)

    @property
    def read_noise_electrons(self) -> float:
        return self.detector.read_noise_electrons_per_pixel_per_readout

    @property
    def camera_shape(self) -> tuple[int, int]:
        return self.grid.camera_shape

    @property
    def roi_pixel_count(self) -> int:
        return self.grid.roi_pixel_count

    def column_density(self, parameters: SmoothTFParameters) -> NDArray[np.floating]:
        return smooth_tf_column_density(self.grid.y_grid_m, self.grid.z_grid_m, parameters)

    def _validate_density_map(self, column_density_m2: ArrayLike) -> NDArray[np.floating]:
        density = np.asarray(column_density_m2, dtype=float)
        if density.shape != self.grid.y_grid_m.shape:
            raise ValueError(
                f"column density must have reconstruction-grid shape {self.grid.y_grid_m.shape}"
            )
        if np.any(~np.isfinite(density)):
            raise ValueError("column density must be finite")
        if np.any(density < 0):
            raise ValueError("column density must be non-negative")
        return density

    def _validate_density_derivatives(
        self,
        density_derivatives_m2: ArrayLike,
    ) -> NDArray[np.floating]:
        derivatives = np.asarray(density_derivatives_m2, dtype=float)
        expected_tail = self.grid.y_grid_m.shape
        if derivatives.ndim != 3 or derivatives.shape[1:] != expected_tail:
            raise ValueError(
                "density derivatives must have shape "
                f"(n_parameter, {expected_tail[0]}, {expected_tail[1]})"
            )
        if np.any(~np.isfinite(derivatives)):
            raise ValueError("column-density derivatives must be finite")
        return derivatives

    def _propagated_linear_fields_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.complexfloating], NDArray[np.complexfloating]]:
        density = self._validate_density_map(column_density_m2)
        theta = self.response.rotation_per_column_density_rad_m2 * density
        sigma_plus = 1.0 + propagate_scattered_field(np.exp(1j * theta) - 1.0, self.grid.pupil)
        sigma_minus = 1.0 + propagate_scattered_field(np.exp(-1j * theta) - 1.0, self.grid.pupil)
        parallel = (sigma_plus + sigma_minus) / 2.0
        perpendicular = 1j * (sigma_plus - sigma_minus) / 2.0
        return parallel, perpendicular

    def _propagated_linear_fields(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.complexfloating], NDArray[np.complexfloating]]:
        return self._propagated_linear_fields_from_density(self.column_density(parameters))

    def _propagated_linear_fields_and_density_jacobian(
        self,
        column_density_m2: ArrayLike,
        density_derivatives_m2: ArrayLike,
    ) -> tuple[
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
    ]:
        parallel, perpendicular = self._propagated_linear_fields_from_density(
            column_density_m2
        )
        derivative_parallel, derivative_perpendicular = self._propagated_derivative_fields(
            column_density_m2,
            density_derivatives_m2,
        )
        return parallel, perpendicular, derivative_parallel, derivative_perpendicular

    def _propagated_derivative_fields(
        self,
        column_density_m2: ArrayLike,
        density_derivatives_m2: ArrayLike,
    ) -> tuple[NDArray[np.complexfloating], NDArray[np.complexfloating]]:
        """Propagate one bounded batch of object-density derivatives."""

        density = self._validate_density_map(column_density_m2)
        density_derivatives = self._validate_density_derivatives(density_derivatives_m2)
        parameter_count = density_derivatives.shape[0]
        coefficient = self.response.rotation_per_column_density_rad_m2
        theta = coefficient * density
        theta_derivatives = coefficient * density_derivatives
        plus_object = np.exp(1j * theta)
        minus_object = np.exp(-1j * theta)
        fields = np.concatenate(
            [
                1j * plus_object[None, ...] * theta_derivatives,
                -1j * minus_object[None, ...] * theta_derivatives,
            ],
            axis=0,
        )
        propagated = scipy_fft.ifft2(
            scipy_fft.fft2(fields, axes=(-2, -1), workers=-1)
            * self.grid.pupil[None, ...],
            axes=(-2, -1),
            workers=-1,
        )
        derivative_plus = propagated[:parameter_count]
        derivative_minus = propagated[parameter_count : 2 * parameter_count]
        derivative_parallel = (derivative_plus + derivative_minus) / 2.0
        derivative_perpendicular = 1j * (derivative_plus - derivative_minus) / 2.0
        return derivative_parallel, derivative_perpendicular

    def _propagated_linear_fields_and_jacobian(
        self,
        internal_vector: ArrayLike,
    ) -> tuple[
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
        NDArray[np.complexfloating],
    ]:
        parameters = from_internal(internal_vector)
        density, density_derivatives = smooth_tf_density_and_internal_jacobian(
            self.grid.y_grid_m,
            self.grid.z_grid_m,
            parameters,
        )
        return self._propagated_linear_fields_and_density_jacobian(
            density,
            density_derivatives,
        )

    def flatten_observed(self, *channels: ArrayLike) -> NDArray[np.floating]:
        if len(channels) != len(self.channel_names):
            raise ValueError(
                f"{type(self).__name__} expects {len(self.channel_names)} observed channels"
            )
        flattened: list[NDArray[np.floating]] = []
        for channel in channels:
            array = np.asarray(channel, dtype=float)
            if array.shape != self.camera_shape:
                raise ValueError(f"observed channel must have shape {self.camera_shape}")
            flattened.append(array[self.grid.roi_mask])
        return np.concatenate(flattened)

    def simulate_channels(
        self,
        parameters: SmoothTFParameters,
        rng: np.random.Generator,
    ) -> tuple[NDArray[np.floating], ...]:
        return self.simulate_channels_from_density(self.column_density(parameters), rng)

    def simulate_channels_from_density(
        self,
        column_density_m2: ArrayLike,
        rng: np.random.Generator,
    ) -> tuple[NDArray[np.floating], ...]:
        return tuple(
            simulate_poisson_gaussian_counts(
                expected,
                read_noise_electrons=self.read_noise_electrons,
                rng=rng,
            )
            for expected in self.expected_channels_from_density(column_density_m2)
        )

    def expected_channels(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.floating], ...]:
        return self.expected_channels_from_density(self.column_density(parameters))

    def expected_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], ...]:
        raise NotImplementedError

    def expected_vector_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> NDArray[np.floating]:
        return self.flatten_observed(*self.expected_channels_from_density(column_density_m2))

    def expected_vector_and_jacobian_from_density(
        self,
        column_density_m2: ArrayLike,
        density_derivatives_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        raise NotImplementedError

    def expected_vector_and_jacobian_from_density_batches(
        self,
        column_density_m2: ArrayLike,
        derivative_batches: Iterable[tuple[slice, NDArray[np.floating]]],
        parameter_count: int,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        raise NotImplementedError

    def expected_vector_and_jacobian_model(
        self,
        model: DifferentiableColumnDensityModel,
        parameter_vector: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        """Evaluate any differentiable object model through the same instrument."""

        density = model.column_density(parameter_vector)
        batches = model.iter_column_density_jacobian(
            parameter_vector,
            self.jacobian_batch_size,
        )
        return self.expected_vector_and_jacobian_from_density_batches(
            density,
            batches,
            model.parameter_count,
        )

    def expected_vector_and_jacobian_internal(
        self,
        internal_vector: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        raise NotImplementedError


class DualPortFaradayMeasurement(_FaradayMeasurementBase):
    """Finite-aperture dual-port model fitted to the raw H and V counts."""

    channel_names = ("H", "V")

    def normalised_channels(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        return self.normalised_channels_from_density(self.column_density(parameters))

    def normalised_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        parallel, perpendicular = self._propagated_linear_fields_from_density(column_density_m2)
        h_field = parallel - perpendicular
        v_field = parallel + perpendicular
        h = self.grid.camera_average(np.abs(h_field) ** 2 / 2.0)
        v = self.grid.camera_average(np.abs(v_field) ** 2 / 2.0)
        return h, v

    def expected_channels(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        return self.expected_channels_from_density(self.column_density(parameters))

    def expected_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        scale = self.detector.photoelectrons_per_i0_pixel
        h, v = self.normalised_channels_from_density(column_density_m2)
        return scale * h, scale * v

    def expected_vector_and_jacobian_internal(
        self,
        internal_vector: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        parameters = from_internal(internal_vector)
        density, derivatives = smooth_tf_density_and_internal_jacobian(
            self.grid.y_grid_m,
            self.grid.z_grid_m,
            parameters,
        )
        return self.expected_vector_and_jacobian_from_density(density, derivatives)

    def expected_vector_and_jacobian_from_density(
        self,
        column_density_m2: ArrayLike,
        density_derivatives_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        derivatives = self._validate_density_derivatives(density_derivatives_m2)

        def batches() -> Iterable[tuple[slice, NDArray[np.floating]]]:
            for start in range(0, derivatives.shape[0], self.jacobian_batch_size):
                stop = min(start + self.jacobian_batch_size, derivatives.shape[0])
                yield slice(start, stop), derivatives[start:stop]

        return self.expected_vector_and_jacobian_from_density_batches(
            column_density_m2,
            batches(),
            derivatives.shape[0],
        )

    def expected_vector_and_jacobian_from_density_batches(
        self,
        column_density_m2: ArrayLike,
        derivative_batches: Iterable[tuple[slice, NDArray[np.floating]]],
        parameter_count: int,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        parallel, perpendicular = self._propagated_linear_fields_from_density(
            column_density_m2
        )
        h_field = parallel - perpendicular
        v_field = parallel + perpendicular
        h = self.grid.camera_average(np.abs(h_field) ** 2 / 2.0)
        v = self.grid.camera_average(np.abs(v_field) ** 2 / 2.0)
        scale = self.detector.photoelectrons_per_i0_pixel
        prediction = np.concatenate(
            [scale * h[self.grid.roi_mask], scale * v[self.grid.roi_mask]]
        )
        jacobian = np.empty((prediction.size, parameter_count), dtype=float)
        populated = np.zeros(parameter_count, dtype=bool)
        for parameter_slice, density_derivatives in derivative_batches:
            if parameter_slice.start is None or parameter_slice.stop is None:
                raise ValueError("Jacobian batch slices must have explicit bounds")
            if (
                parameter_slice.start < 0
                or parameter_slice.stop > parameter_count
                or parameter_slice.start >= parameter_slice.stop
            ):
                raise ValueError("Jacobian batch slice lies outside the parameter vector")
            if np.any(populated[parameter_slice]):
                raise ValueError("Jacobian batches overlap")
            d_parallel, d_perpendicular = self._propagated_derivative_fields(
                column_density_m2,
                density_derivatives,
            )
            expected_batch = parameter_slice.stop - parameter_slice.start
            if d_parallel.shape[0] != expected_batch:
                raise ValueError("Jacobian batch width does not match its parameter slice")
            d_h_field = d_parallel - d_perpendicular
            d_v_field = d_parallel + d_perpendicular
            d_h = self.grid.camera_average_stack(
                np.real(np.conj(h_field)[None, ...] * d_h_field)
            )
            d_v = self.grid.camera_average_stack(
                np.real(np.conj(v_field)[None, ...] * d_v_field)
            )
            jacobian[:, parameter_slice] = np.concatenate(
                [scale * d_h[:, self.grid.roi_mask], scale * d_v[:, self.grid.roi_mask]],
                axis=1,
            ).T
            populated[parameter_slice] = True
        if not np.all(populated):
            raise ValueError("Jacobian batches do not cover every parameter")
        return prediction, jacobian


class DarkFieldFaradayMeasurement(_FaradayMeasurementBase):
    """Finite-aperture crossed-analyser model fitted to one raw count channel."""

    channel_names = ("dark_field",)

    def normalised_channels(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.floating]]:
        return self.normalised_channels_from_density(self.column_density(parameters))

    def normalised_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating]]:
        _, perpendicular = self._propagated_linear_fields_from_density(column_density_m2)
        image = self.grid.camera_average(np.abs(perpendicular) ** 2)
        return (image,)

    def expected_channels(
        self,
        parameters: SmoothTFParameters,
    ) -> tuple[NDArray[np.floating]]:
        return self.expected_channels_from_density(self.column_density(parameters))

    def expected_channels_from_density(
        self,
        column_density_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating]]:
        scale = self.detector.photoelectrons_per_i0_pixel
        return tuple(
            scale * image for image in self.normalised_channels_from_density(column_density_m2)
        )

    def expected_vector_and_jacobian_internal(
        self,
        internal_vector: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        parameters = from_internal(internal_vector)
        density, derivatives = smooth_tf_density_and_internal_jacobian(
            self.grid.y_grid_m,
            self.grid.z_grid_m,
            parameters,
        )
        return self.expected_vector_and_jacobian_from_density(density, derivatives)

    def expected_vector_and_jacobian_from_density(
        self,
        column_density_m2: ArrayLike,
        density_derivatives_m2: ArrayLike,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        derivatives = self._validate_density_derivatives(density_derivatives_m2)

        def batches() -> Iterable[tuple[slice, NDArray[np.floating]]]:
            for start in range(0, derivatives.shape[0], self.jacobian_batch_size):
                stop = min(start + self.jacobian_batch_size, derivatives.shape[0])
                yield slice(start, stop), derivatives[start:stop]

        return self.expected_vector_and_jacobian_from_density_batches(
            column_density_m2,
            batches(),
            derivatives.shape[0],
        )

    def expected_vector_and_jacobian_from_density_batches(
        self,
        column_density_m2: ArrayLike,
        derivative_batches: Iterable[tuple[slice, NDArray[np.floating]]],
        parameter_count: int,
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        _, perpendicular = self._propagated_linear_fields_from_density(column_density_m2)
        image = self.grid.camera_average(np.abs(perpendicular) ** 2)
        scale = self.detector.photoelectrons_per_i0_pixel
        prediction = scale * image[self.grid.roi_mask]
        jacobian = np.empty((prediction.size, parameter_count), dtype=float)
        populated = np.zeros(parameter_count, dtype=bool)
        for parameter_slice, density_derivatives in derivative_batches:
            if parameter_slice.start is None or parameter_slice.stop is None:
                raise ValueError("Jacobian batch slices must have explicit bounds")
            if (
                parameter_slice.start < 0
                or parameter_slice.stop > parameter_count
                or parameter_slice.start >= parameter_slice.stop
            ):
                raise ValueError("Jacobian batch slice lies outside the parameter vector")
            if np.any(populated[parameter_slice]):
                raise ValueError("Jacobian batches overlap")
            _, d_perpendicular = self._propagated_derivative_fields(
                column_density_m2,
                density_derivatives,
            )
            expected_batch = parameter_slice.stop - parameter_slice.start
            if d_perpendicular.shape[0] != expected_batch:
                raise ValueError("Jacobian batch width does not match its parameter slice")
            derivatives = self.grid.camera_average_stack(
                2.0 * np.real(np.conj(perpendicular)[None, ...] * d_perpendicular)
            )
            jacobian[:, parameter_slice] = (
                scale * derivatives[:, self.grid.roi_mask]
            ).T
            populated[parameter_slice] = True
        if not np.all(populated):
            raise ValueError("Jacobian batches do not cover every parameter")
        return prediction, jacobian


def assert_response_scale_degeneracy(
    measurement: _FaradayMeasurementBase,
    parameters: SmoothTFParameters,
    *,
    scale: float,
    atol: float = 1e-10,
) -> None:
    """Check the exact synthetic degeneracy between density scale and kappa_F.

    This helper is intended for tests and documentation. It raises if replacing
    ``kappa_F`` by ``scale*kappa_F`` and the density amplitude by its inverse
    changes the predicted count channels.
    """

    if scale <= 0:
        raise ValueError("scale must be positive")
    scaled_response = FaradayResponseContract(
        phase_per_column_density_rad_m2=measurement.response.phase_per_column_density_rad_m2,
        kappa_f=measurement.response.kappa_f * scale,
    )
    scaled_parameters = SmoothTFParameters(
        column_density_peak_m2=parameters.column_density_peak_m2 / scale,
        y0_um=parameters.y0_um,
        z0_um=parameters.z0_um,
        radius_y_um=parameters.radius_y_um,
        radius_z_um=parameters.radius_z_um,
    )
    scaled_measurement = type(measurement)(
        grid=measurement.grid,
        detector=measurement.detector,
        response=scaled_response,
    )
    for original, transformed in zip(
        measurement.expected_channels(parameters),
        scaled_measurement.expected_channels(scaled_parameters),
        strict=True,
    ):
        np.testing.assert_allclose(original, transformed, rtol=0.0, atol=atol)
