"""Source-declared analytic truth maps for the initial-condition suite.

The generators in this module translate the parameters declared in
``dpfi_initial_condition_suite_v1_orca_fusion_m10.json`` into two-dimensional
column-density maps in ``m^-2``.  They are deliberately separate from the
inverse basis, optical forward operator and detector-noise model.

The coordinate convention is the one used by the active reconstruction: the
probe propagates along the first trap axis (``x``), while the supplied object
grids span the second and third axes (``y`` and ``z``).  Contact Thomas--Fermi
profiles are analytic reference/surrogate states.  The ideal-gas halo is a
population-normalised saturated Bose profile, and the modulated and Gaussian
maps are declared density surrogates rather than equilibrium dipolar states.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import special

from ..atomic_model import ThomasFermiState, build_thomas_fermi_state
from ..profiles import thomas_fermi_profile_2d
from .synthetic_morphologies import SyntheticMorphology


MetadataValue: TypeAlias = str | int | float | bool | tuple[float, ...]


@dataclass(frozen=True)
class InitialConditionTruth:
    """One immutable initial-condition truth map and its provenance.

    ``component_items`` records the physical/analytic components before they
    are summed.  It is intended for provenance and observable interpretation;
    reconstruction consumes only ``morphology``.
    """

    condition_id: str
    generator_name: str
    morphology: SyntheticMorphology
    component_items: tuple[tuple[str, NDArray[np.floating]], ...]
    source_items: tuple[tuple[str, MetadataValue], ...]
    derived_items: tuple[tuple[str, MetadataValue], ...]

    def __post_init__(self) -> None:
        if not self.condition_id:
            raise ValueError("initial-condition id must be non-empty")
        if not self.generator_name:
            raise ValueError("initial-condition generator name must be non-empty")
        if self.morphology.name != self.condition_id:
            raise ValueError("morphology name must equal the initial-condition id")

        total = _immutable_density(
            self.morphology.column_density_m2,
            name="initial-condition density",
            require_nonempty=True,
        )
        morphology = SyntheticMorphology(
            name=self.morphology.name,
            column_density_m2=total,
            description=self.morphology.description,
            feature_class=self.morphology.feature_class,
        )

        component_items: list[tuple[str, NDArray[np.floating]]] = []
        component_names: list[str] = []
        component_sum = np.zeros_like(total)
        for component_name, values in self.component_items:
            if not component_name:
                raise ValueError("initial-condition component names must be non-empty")
            component = _immutable_density(
                values,
                name=f"initial-condition component {component_name!r}",
                require_nonempty=False,
            )
            if component.shape != total.shape:
                raise ValueError(
                    "initial-condition components must match the total density shape"
                )
            component_names.append(component_name)
            component_items.append((component_name, component))
            component_sum += component
        if not component_items or len(component_names) != len(set(component_names)):
            raise ValueError(
                "initial-condition components must be non-empty and uniquely named"
            )
        if not np.array_equal(component_sum, total):
            raise ValueError("initial-condition components must sum exactly to the total")

        source_items = tuple(self.source_items)
        derived_items = tuple(self.derived_items)
        _validate_metadata_items(source_items, "source")
        _validate_metadata_items(derived_items, "derived")
        object.__setattr__(self, "morphology", morphology)
        object.__setattr__(self, "component_items", tuple(component_items))
        object.__setattr__(self, "source_items", source_items)
        object.__setattr__(self, "derived_items", derived_items)

    @property
    def column_density_m2(self) -> NDArray[np.floating]:
        """Return the immutable total column-density map in ``m^-2``."""

        return self.morphology.column_density_m2

    @property
    def description(self) -> str:
        return self.morphology.description

    @property
    def feature_class(self) -> str:
        return self.morphology.feature_class

    @property
    def components(self) -> dict[str, NDArray[np.floating]]:
        """Return a new mapping to the immutable component arrays."""

        return dict(self.component_items)

    @property
    def source_metadata(self) -> dict[str, MetadataValue]:
        return dict(self.source_items)

    @property
    def derived_metadata(self) -> dict[str, MetadataValue]:
        return dict(self.derived_items)


@dataclass(frozen=True)
class _AtomicInputs:
    mass_kg: float
    bohr_radius_m: float
    hbar_j_s: float
    boltzmann_constant_j_k: float


@dataclass(frozen=True)
class _ConditionIdentity:
    condition_id: str
    generator_name: str
    description: str
    feature_class: str


_SOURCE_KEYS = (
    "family",
    "parameter_status",
    "source_label",
    "source_url",
    "source_data_url",
    "relationship_to_existing_reference",
)


def _immutable_density(
    values: ArrayLike,
    *,
    name: str,
    require_nonempty: bool,
) -> NDArray[np.floating]:
    density = np.asarray(values, dtype=float)
    if density.ndim != 2 or np.any(~np.isfinite(density)):
        raise ValueError(f"{name} must be a finite two-dimensional map")
    if np.any(density < 0.0):
        raise ValueError(f"{name} must be non-negative")
    if require_nonempty and not np.any(density > 0.0):
        raise ValueError(f"{name} must be non-empty")
    result = np.array(density, copy=True)
    result.setflags(write=False)
    return result


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
        values = value if isinstance(value, tuple) else (value,)
        if any(not np.isfinite(float(item)) for item in values):
            raise ValueError(f"{kind} metadata value for {key!r} must be finite")


def _string(mapping: Mapping[str, object], key: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value:
        raise ValueError(f"{key} must be a non-empty string")
    return value


def _number(
    mapping: Mapping[str, object],
    key: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
) -> float:
    value = mapping.get(key)
    if isinstance(value, bool):
        raise ValueError(f"{key} must be numeric")
    try:
        checked = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key} must be numeric") from exc
    if not np.isfinite(checked):
        raise ValueError(f"{key} must be finite")
    if positive and checked <= 0.0:
        raise ValueError(f"{key} must be positive")
    if nonnegative and checked < 0.0:
        raise ValueError(f"{key} must be non-negative")
    return checked


def _number_tuple(
    mapping: Mapping[str, object],
    key: str,
    *,
    length: int,
    positive: bool = False,
    nonnegative: bool = False,
) -> tuple[float, ...]:
    raw = mapping.get(key)
    if isinstance(raw, (str, bytes)):
        raise ValueError(f"{key} must contain {length} numeric values")
    try:
        raw_values = tuple(raw)  # type: ignore[arg-type]
    except TypeError as exc:
        raise ValueError(f"{key} must contain {length} numeric values") from exc
    if len(raw_values) != length:
        raise ValueError(f"{key} must contain {length} numeric values")
    values: list[float] = []
    for value in raw_values:
        if isinstance(value, bool):
            raise ValueError(f"{key} values must be numeric")
        try:
            checked = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} values must be numeric") from exc
        if not np.isfinite(checked):
            raise ValueError(f"{key} values must be finite")
        if positive and checked <= 0.0:
            raise ValueError(f"{key} values must be positive")
        if nonnegative and checked < 0.0:
            raise ValueError(f"{key} values must be non-negative")
        values.append(checked)
    return tuple(values)


def _validated_grids(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
) -> tuple[NDArray[np.floating], NDArray[np.floating], float]:
    y = np.asarray(y_grid_m, dtype=float)
    z = np.asarray(z_grid_m, dtype=float)
    if y.ndim != 2 or y.shape != z.shape:
        raise ValueError("initial-condition grids must be same-shape 2D arrays")
    if min(y.shape) < 2 or np.any(~np.isfinite(y)) or np.any(~np.isfinite(z)):
        raise ValueError("initial-condition grids must be finite and at least 2 by 2")
    y_axis = y[0, :]
    z_axis = z[:, 0]
    scale = max(float(np.max(np.abs(y))), float(np.max(np.abs(z))), 1.0)
    tolerance = 32.0 * np.finfo(float).eps * scale
    if not np.allclose(y, y_axis[None, :], rtol=0.0, atol=tolerance):
        raise ValueError("initial-condition y grid must be a rectilinear meshgrid")
    if not np.allclose(z, z_axis[:, None], rtol=0.0, atol=tolerance):
        raise ValueError("initial-condition z grid must be a rectilinear meshgrid")
    dy = np.diff(y_axis)
    dz = np.diff(z_axis)
    if np.any(dy <= 0.0) or np.any(dz <= 0.0):
        raise ValueError("initial-condition coordinate axes must be strictly increasing")
    if not np.allclose(dy, dy[0], rtol=1e-10, atol=0.0) or not np.allclose(
        dz, dz[0], rtol=1e-10, atol=0.0
    ):
        raise ValueError("initial-condition coordinate axes must be uniformly spaced")
    return y, z, float(dy[0] * dz[0])


def _atomic_inputs(atomic_constants: Mapping[str, object]) -> _AtomicInputs:
    mass_number = _number(atomic_constants, "mass_number", positive=True)
    atomic_mass_unit = _number(
        atomic_constants, "atomic_mass_unit_kg", positive=True
    )
    return _AtomicInputs(
        mass_kg=mass_number * atomic_mass_unit,
        bohr_radius_m=_number(atomic_constants, "bohr_radius_m", positive=True),
        hbar_j_s=_number(atomic_constants, "hbar_j_s", positive=True),
        boltzmann_constant_j_k=_number(
            atomic_constants, "boltzmann_constant_j_k", positive=True
        ),
    )


def _identity(
    condition: Mapping[str, object], expected_generator: str
) -> _ConditionIdentity:
    generator = _string(condition, "generator")
    if generator != expected_generator:
        raise ValueError(
            f"expected generator {expected_generator!r}, received {generator!r}"
        )
    return _ConditionIdentity(
        condition_id=_string(condition, "id"),
        generator_name=generator,
        description=_string(condition, "description"),
        feature_class=_string(condition, "feature_class"),
    )


def _source_items(condition: Mapping[str, object]) -> tuple[tuple[str, MetadataValue], ...]:
    items: list[tuple[str, MetadataValue]] = []
    for key in _SOURCE_KEYS:
        if key in condition:
            value = condition[key]
            if not isinstance(value, str):
                raise ValueError(f"{key} must be a string when supplied")
            items.append((key, value))
    return tuple(items)


def _tf_component(
    y_grid_m: NDArray[np.floating],
    z_grid_m: NDArray[np.floating],
    *,
    atom_number: float,
    scattering_length_bohr: float,
    trap_frequencies_hz: tuple[float, float, float],
    atomic: _AtomicInputs,
) -> tuple[NDArray[np.floating], ThomasFermiState]:
    if atom_number <= 0.0:
        raise ValueError("contact-TF atom number must be positive")
    state = build_thomas_fermi_state(
        atom_number,
        scattering_length_bohr * atomic.bohr_radius_m,
        trap_frequencies_hz,
        atomic.mass_kg,
        atomic.hbar_j_s,
        atomic.boltzmann_constant_j_k,
    )
    density = state.column_density[0] * thomas_fermi_profile_2d(
        y_grid_m,
        z_grid_m,
        state.radii[1],
        state.radii[2],
    )
    return np.asarray(density, dtype=float), state


def _tf_metadata(state: ThomasFermiState) -> dict[str, MetadataValue]:
    return {
        "tf_radii_um": tuple(float(value * 1e6) for value in state.radii),
        "tf_peak_3d_density_m3": float(state.peak_density),
        "tf_peak_column_density_m2": float(state.column_density[0]),
        "tf_chemical_potential_j": float(state.chemical_potential),
        "tf_chemical_potential_temperature_nk": float(
            state.chemical_potential_temperature * 1e9
        ),
        "tf_atom_number_check": float(state.atom_number_check),
    }


def _zero_tf_metadata() -> dict[str, MetadataValue]:
    return {
        "tf_radii_um": (0.0, 0.0, 0.0),
        "tf_peak_3d_density_m3": 0.0,
        "tf_peak_column_density_m2": 0.0,
        "tf_chemical_potential_j": 0.0,
        "tf_chemical_potential_temperature_nk": 0.0,
        "tf_atom_number_check": 0.0,
    }


def _ideal_bose_halo(
    y_grid_m: NDArray[np.floating],
    z_grid_m: NDArray[np.floating],
    *,
    atom_number: float,
    temperature_nk: float,
    trap_frequencies_hz: tuple[float, float, float],
    atomic: _AtomicInputs,
) -> tuple[NDArray[np.floating], dict[str, MetadataValue]]:
    if atom_number == 0.0:
        return np.zeros_like(y_grid_m), {
            "ideal_bose_thermal_scale_y_um": 0.0,
            "ideal_bose_thermal_scale_z_um": 0.0,
            "ideal_bose_critical_atoms": 0.0,
            "ideal_bose_population_to_critical_ratio": 0.0,
        }
    if temperature_nk <= 0.0:
        raise ValueError("temperature_nk must be positive when thermal_atoms is non-zero")

    temperature_k = temperature_nk * 1e-9
    omega_x, omega_y, omega_z = 2.0 * np.pi * np.asarray(
        trap_frequencies_hz, dtype=float
    )
    beta = 1.0 / (atomic.boltzmann_constant_j_k * temperature_k)
    transverse_potential = 0.5 * atomic.mass_kg * (
        omega_y**2 * y_grid_m**2 + omega_z**2 * z_grid_m**2
    )
    bose_argument = np.exp(-beta * transverse_potential)
    # scipy.special.spence(q) = Li_2(1-q), hence this is Li_2(exp(-beta V)).
    bose_g2 = special.spence(1.0 - bose_argument)
    zeta_three = float(special.zeta(3.0, 1.0))
    prefactor_m2 = (
        atom_number
        * beta
        * atomic.mass_kg
        * omega_y
        * omega_z
        / (2.0 * np.pi * zeta_three)
    )
    density = prefactor_m2 * bose_g2
    omega_bar = float((omega_x * omega_y * omega_z) ** (1.0 / 3.0))
    critical_atoms = zeta_three * (
        atomic.boltzmann_constant_j_k
        * temperature_k
        / (atomic.hbar_j_s * omega_bar)
    ) ** 3
    return np.asarray(density, dtype=float), {
        "ideal_bose_thermal_scale_y_um": float(
            np.sqrt(1.0 / (beta * atomic.mass_kg * omega_y**2)) * 1e6
        ),
        "ideal_bose_thermal_scale_z_um": float(
            np.sqrt(1.0 / (beta * atomic.mass_kg * omega_z**2)) * 1e6
        ),
        "ideal_bose_critical_atoms": float(critical_atoms),
        "ideal_bose_population_to_critical_ratio": float(
            atom_number / critical_atoms
        ),
    }


def _common_inputs(
    condition: Mapping[str, object],
) -> tuple[
    float,
    float,
    float,
    tuple[float, float, float],
    float,
    float,
]:
    condensate_atoms = _number(condition, "condensate_atoms", nonnegative=True)
    thermal_atoms = _number(condition, "thermal_atoms", nonnegative=True)
    temperature_nk = _number(condition, "temperature_nk", nonnegative=True)
    trap = _number_tuple(
        condition, "trap_frequencies_hz", length=3, positive=True
    )
    scattering_length_bohr = _number(
        condition, "scattering_length_bohr", positive=True
    )
    dipole_orientation_deg = _number(condition, "dipole_orientation_deg")
    return (
        condensate_atoms,
        thermal_atoms,
        temperature_nk,
        (trap[0], trap[1], trap[2]),
        scattering_length_bohr,
        dipole_orientation_deg,
    )


def _truth(
    *,
    identity: _ConditionIdentity,
    condition: Mapping[str, object],
    components: tuple[tuple[str, NDArray[np.floating]], ...],
    cell_area_m2: float,
    declared_condensate_atoms: float,
    declared_thermal_atoms: float,
    derived: Mapping[str, MetadataValue],
) -> InitialConditionTruth:
    total = np.zeros_like(components[0][1], dtype=float)
    captured_items: list[tuple[str, MetadataValue]] = []
    for component_name, component in components:
        total += component
        captured_items.append(
            (
                f"grid_captured_{component_name}_atoms",
                float(np.sum(component) * cell_area_m2),
            )
        )
    total.setflags(write=False)
    declared_total = declared_condensate_atoms + declared_thermal_atoms
    captured_total = float(np.sum(total) * cell_area_m2)
    derived_items: list[tuple[str, MetadataValue]] = [
        ("line_of_sight_axis", "x"),
        ("object_plane_axes", "y,z"),
        ("declared_condensate_atoms", declared_condensate_atoms),
        ("declared_thermal_atoms", declared_thermal_atoms),
        ("declared_total_atoms", declared_total),
        ("grid_cell_area_m2", cell_area_m2),
        *captured_items,
        ("grid_captured_total_atoms", captured_total),
        ("grid_captured_fraction_of_declared_total", captured_total / declared_total),
        ("grid_peak_column_density_m2", float(np.max(total))),
        *derived.items(),
    ]
    return InitialConditionTruth(
        condition_id=identity.condition_id,
        generator_name=identity.generator_name,
        morphology=SyntheticMorphology(
            identity.condition_id,
            total,
            identity.description,
            identity.feature_class,
        ),
        component_items=components,
        source_items=_source_items(condition),
        derived_items=tuple(derived_items),
    )


def build_contact_tf_truth(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    condition: Mapping[str, object],
    atomic_constants: Mapping[str, object],
) -> InitialConditionTruth:
    """Build a contact-interaction Thomas--Fermi column-density truth."""

    identity = _identity(condition, "contact_tf")
    y_grid, z_grid, cell_area = _validated_grids(y_grid_m, z_grid_m)
    atomic = _atomic_inputs(atomic_constants)
    n0, nth, temperature_nk, trap, scattering_bohr, dipole_angle = _common_inputs(
        condition
    )
    if n0 <= 0.0:
        raise ValueError("contact_tf requires a positive condensate_atoms value")
    if nth != 0.0:
        raise ValueError("contact_tf does not represent a non-zero thermal component")
    condensate, state = _tf_component(
        y_grid,
        z_grid,
        atom_number=n0,
        scattering_length_bohr=scattering_bohr,
        trap_frequencies_hz=trap,
        atomic=atomic,
    )
    derived: dict[str, MetadataValue] = {
        "trap_frequencies_hz": trap,
        "scattering_length_bohr": scattering_bohr,
        "temperature_nk": temperature_nk,
        "dipole_orientation_deg": dipole_angle,
        "density_model": "contact_thomas_fermi",
        "dipolar_mean_field_included": False,
        **_tf_metadata(state),
    }
    return _truth(
        identity=identity,
        condition=condition,
        components=(("condensate_core", condensate),),
        cell_area_m2=cell_area,
        declared_condensate_atoms=n0,
        declared_thermal_atoms=nth,
        derived=derived,
    )


def build_bimodal_ideal_gas_truth(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    condition: Mapping[str, object],
    atomic_constants: Mapping[str, object],
) -> InitialConditionTruth:
    """Build a contact-TF core plus a population-normalised ideal-Bose halo.

    ``condensate_atoms`` and ``thermal_atoms`` are direct inputs.  The fixed
    fugacity controls only the analytic halo shape; it is not used to infer a
    total population or impose a thermal-saturation closure.
    """

    identity = _identity(condition, "bimodal_ideal_gas")
    y_grid, z_grid, cell_area = _validated_grids(y_grid_m, z_grid_m)
    atomic = _atomic_inputs(atomic_constants)
    n0, nth, temperature_nk, trap, scattering_bohr, dipole_angle = _common_inputs(
        condition
    )
    if n0 + nth <= 0.0:
        raise ValueError("bimodal_ideal_gas requires a non-zero declared population")
    if n0 > 0.0:
        condensate, state = _tf_component(
            y_grid,
            z_grid,
            atom_number=n0,
            scattering_length_bohr=scattering_bohr,
            trap_frequencies_hz=trap,
            atomic=atomic,
        )
        tf_metadata = _tf_metadata(state)
    else:
        condensate = np.zeros_like(y_grid)
        tf_metadata = _zero_tf_metadata()
    thermal, thermal_metadata = _ideal_bose_halo(
        y_grid,
        z_grid,
        atom_number=nth,
        temperature_nk=temperature_nk,
        trap_frequencies_hz=trap,
        atomic=atomic,
    )
    derived: dict[str, MetadataValue] = {
        "trap_frequencies_hz": trap,
        "scattering_length_bohr": scattering_bohr,
        "temperature_nk": temperature_nk,
        "dipole_orientation_deg": dipole_angle,
        "density_model": "contact_tf_plus_population_normalised_ideal_bose_shape",
        "thermal_fugacity": 1.0,
        "dipolar_mean_field_included": False,
        **tf_metadata,
        **thermal_metadata,
    }
    return _truth(
        identity=identity,
        condition=condition,
        components=(
            ("condensate_core", condensate),
            ("thermal_halo", thermal),
        ),
        cell_area_m2=cell_area,
        declared_condensate_atoms=n0,
        declared_thermal_atoms=nth,
        derived=derived,
    )


def _normalised_bessel_transform(order: float, argument: float) -> float:
    """Return Gamma(order+1) (2/x)^order J_order(x), stable near zero."""

    if abs(argument) < 1e-4:
        return float(
            1.0
            - argument**2 / (4.0 * (order + 1.0))
            + argument**4 / (32.0 * (order + 1.0) * (order + 2.0))
        )
    return float(
        special.gamma(order + 1.0)
        * (2.0 / argument) ** order
        * special.jv(order, argument)
    )


def build_tf_roton_modulation_truth(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    condition: Mapping[str, object],
    atomic_constants: Mapping[str, object],
) -> InitialConditionTruth:
    """Build a number-preserving modulated contact-TF core plus Bose halo."""

    identity = _identity(condition, "tf_roton_modulation")
    y_grid, z_grid, cell_area = _validated_grids(y_grid_m, z_grid_m)
    atomic = _atomic_inputs(atomic_constants)
    n0, nth, temperature_nk, trap, scattering_bohr, dipole_angle = _common_inputs(
        condition
    )
    if n0 <= 0.0:
        raise ValueError("tf_roton_modulation requires positive condensate_atoms")
    contrast = _number(condition, "modulation_contrast", nonnegative=True)
    if contrast > 1.0:
        raise ValueError("modulation_contrast must lie in [0, 1]")
    wavevector_factor = _number(
        condition, "roton_wavevector_factor", positive=True
    )
    phase = (
        _number(condition, "modulation_phase_rad")
        if "modulation_phase_rad" in condition
        else 0.0
    )
    base_condensate, state = _tf_component(
        y_grid,
        z_grid,
        atom_number=n0,
        scattering_length_bohr=scattering_bohr,
        trap_frequencies_hz=trap,
        atomic=atomic,
    )
    omega_z = 2.0 * np.pi * trap[2]
    oscillator_length_z_m = np.sqrt(
        atomic.hbar_j_s / (atomic.mass_kg * omega_z)
    )
    wavevector_m_inv = wavevector_factor / oscillator_length_z_m
    transform_argument = float(wavevector_m_inv * state.radii[1])
    cosine_mean = np.cos(phase) * _normalised_bessel_transform(
        2.5, transform_argument
    )
    population_normalisation = float(1.0 + contrast * cosine_mean)
    if population_normalisation <= 0.0:
        raise ValueError("TF modulation population normalisation must be positive")
    modulation = 1.0 + contrast * np.cos(wavevector_m_inv * y_grid + phase)
    modulated_condensate = base_condensate * modulation / population_normalisation
    # The high-frequency surrogate is represented only on the declared object
    # grid.  Preserve its configured population under the same cell-area rule
    # used by downstream observables, avoiding a grid-sampling atom-number
    # drift while leaving the analytic modulation and its phase unchanged.
    grid_population_before_correction = float(
        np.sum(modulated_condensate) * cell_area
    )
    grid_population_correction = n0 / grid_population_before_correction
    modulated_condensate *= grid_population_correction
    thermal, thermal_metadata = _ideal_bose_halo(
        y_grid,
        z_grid,
        atom_number=nth,
        temperature_nk=temperature_nk,
        trap_frequencies_hz=trap,
        atomic=atomic,
    )
    derived: dict[str, MetadataValue] = {
        "trap_frequencies_hz": trap,
        "scattering_length_bohr": scattering_bohr,
        "temperature_nk": temperature_nk,
        "dipole_orientation_deg": dipole_angle,
        "density_model": "number_preserving_longitudinal_tf_modulation_plus_ideal_bose",
        "modulation_axis": "y",
        "modulation_contrast": contrast,
        "modulation_phase_rad": phase,
        "roton_wavevector_factor": wavevector_factor,
        "oscillator_length_z_um": float(oscillator_length_z_m * 1e6),
        "modulation_wavevector_um_inv": float(wavevector_m_inv * 1e-6),
        "modulation_period_um": float(2.0 * np.pi / wavevector_m_inv * 1e6),
        "modulation_population_normalisation": population_normalisation,
        "modulation_grid_population_before_correction": (
            grid_population_before_correction
        ),
        "modulation_grid_population_correction": grid_population_correction,
        "dipolar_mean_field_included": False,
        **_tf_metadata(state),
        **thermal_metadata,
    }
    return _truth(
        identity=identity,
        condition=condition,
        components=(
            ("modulated_condensate_core", modulated_condensate),
            ("thermal_halo", thermal),
        ),
        cell_area_m2=cell_area,
        declared_condensate_atoms=n0,
        declared_thermal_atoms=nth,
        derived=derived,
    )


def build_gaussian_chain_truth(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    condition: Mapping[str, object],
    atomic_constants: Mapping[str, object],
) -> InitialConditionTruth:
    """Build an explicitly positioned, continuously population-normalised chain."""

    identity = _identity(condition, "gaussian_chain")
    y_grid, z_grid, cell_area = _validated_grids(y_grid_m, z_grid_m)
    _atomic_inputs(atomic_constants)  # Validate the shared suite contract.
    n0, nth, temperature_nk, trap, scattering_bohr, dipole_angle = _common_inputs(
        condition
    )
    if n0 <= 0.0:
        raise ValueError("gaussian_chain requires positive condensate_atoms")
    if nth != 0.0:
        raise ValueError("gaussian_chain does not represent a thermal component")
    raw_displacements = condition.get("peak_displacements_um")
    raw_weights = condition.get("peak_weights")
    if isinstance(raw_displacements, (str, bytes)) or isinstance(
        raw_weights, (str, bytes)
    ):
        raise ValueError("Gaussian-chain positions and weights must be numeric lists")
    try:
        displacement_values = tuple(raw_displacements)  # type: ignore[arg-type]
        weight_values = tuple(raw_weights)  # type: ignore[arg-type]
    except TypeError as exc:
        raise ValueError(
            "Gaussian-chain positions and weights must be numeric lists"
        ) from exc
    if not displacement_values or len(displacement_values) != len(weight_values):
        raise ValueError(
            "peak_displacements_um and peak_weights must have equal non-zero length"
        )
    displacements = tuple(float(value) for value in displacement_values)
    weights = tuple(float(value) for value in weight_values)
    if (
        any(not np.isfinite(value) for value in displacements)
        or any(not np.isfinite(value) or value < 0.0 for value in weights)
        or not any(value > 0.0 for value in weights)
    ):
        raise ValueError(
            "Gaussian-chain positions must be finite and weights finite/non-negative"
        )
    sigma_y_um = _number(condition, "sigma_y_um", positive=True)
    sigma_z_um = _number(condition, "sigma_z_um", positive=True)
    sigma_y_m = sigma_y_um * 1e-6
    sigma_z_m = sigma_z_um * 1e-6
    values = np.zeros_like(y_grid)
    for displacement_um, weight in zip(displacements, weights, strict=True):
        values += weight * np.exp(
            -0.5
            * (
                ((y_grid - displacement_um * 1e-6) / sigma_y_m) ** 2
                + (z_grid / sigma_z_m) ** 2
            )
        )
    normalisation = n0 / (
        2.0 * np.pi * sigma_y_m * sigma_z_m * sum(weights)
    )
    condensate = normalisation * values
    weighted_centre_um = float(
        np.dot(np.asarray(displacements), np.asarray(weights)) / sum(weights)
    )
    derived: dict[str, MetadataValue] = {
        "trap_frequencies_hz": trap,
        "scattering_length_bohr": scattering_bohr,
        "temperature_nk": temperature_nk,
        "dipole_orientation_deg": dipole_angle,
        "density_model": "population_normalised_gaussian_chain",
        "peak_count": len(displacements),
        "peak_displacements_um": displacements,
        "peak_weights": weights,
        "sigma_y_um": sigma_y_um,
        "sigma_z_um": sigma_z_um,
        "weighted_centre_y_um": weighted_centre_um,
        "dipolar_mean_field_included": False,
        "phase_coherence_included": False,
    }
    return _truth(
        identity=identity,
        condition=condition,
        components=(("condensate_chain", condensate),),
        cell_area_m2=cell_area,
        declared_condensate_atoms=n0,
        declared_thermal_atoms=nth,
        derived=derived,
    )


_GENERATORS = {
    "contact_tf": build_contact_tf_truth,
    "bimodal_ideal_gas": build_bimodal_ideal_gas_truth,
    "tf_roton_modulation": build_tf_roton_modulation_truth,
    "gaussian_chain": build_gaussian_chain_truth,
}


def build_initial_condition_truth(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    condition: Mapping[str, object],
    atomic_constants: Mapping[str, object],
) -> InitialConditionTruth:
    """Dispatch one config condition to its declared analytic truth generator."""

    generator_name = _string(condition, "generator")
    try:
        generator = _GENERATORS[generator_name]
    except KeyError as exc:
        supported = ", ".join(sorted(_GENERATORS))
        raise ValueError(
            f"unsupported initial-condition generator {generator_name!r}; "
            f"expected one of {supported}"
        ) from exc
    return generator(y_grid_m, z_grid_m, condition, atomic_constants)


__all__ = [
    "InitialConditionTruth",
    "build_bimodal_ideal_gas_truth",
    "build_contact_tf_truth",
    "build_gaussian_chain_truth",
    "build_initial_condition_truth",
    "build_tf_roton_modulation_truth",
]
