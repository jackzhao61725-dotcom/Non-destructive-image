"""Tests for independent analytic reconstruction truth inputs."""

from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.reconstruction.morphology_providers import (
    AnalyticTruthSequence,
    build_analytic_truth_catalogue,
    central_hole_cloud,
    controlled_unknown_sequence,
    core_halo_cloud,
    curved_notched_cloud,
    expanded_shifted_cloud,
    gaussian_peak_chain,
    latent_phase_impossibility_pair,
    rotated_asymmetric_cloud,
    transverse_modulated_cloud,
)


PEAK_DENSITY_M2 = 5.0e14


@pytest.fixture(scope="module")
def truth_grid() -> tuple[np.ndarray, np.ndarray]:
    y_axis_m = np.linspace(-40.0e-6, 40.0e-6, 321)
    z_axis_m = np.linspace(-15.0e-6, 15.0e-6, 181)
    return np.meshgrid(y_axis_m, z_axis_m)


def _density_moments(
    density: np.ndarray,
    y_grid_m: np.ndarray,
    z_grid_m: np.ndarray,
) -> tuple[float, float, float, float, float]:
    total = float(np.sum(density))
    centre_y = float(np.sum(density * y_grid_m) / total)
    centre_z = float(np.sum(density * z_grid_m) / total)
    dy = y_grid_m - centre_y
    dz = z_grid_m - centre_z
    return (
        centre_y,
        centre_z,
        float(np.sum(density * dy**2) / total),
        float(np.sum(density * dz**2) / total),
        float(np.sum(density * dy * dz) / total),
    )


def _local_maximum_count(profile: np.ndarray, threshold: float = 0.25) -> int:
    interior = profile[1:-1]
    maxima = (
        (interior > profile[:-2])
        & (interior >= profile[2:])
        & (interior > threshold * np.max(profile))
    )
    return int(np.count_nonzero(maxima))


def test_catalogue_maps_are_nonnegative_peak_normalised_and_distinct(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    catalogue = build_analytic_truth_catalogue(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
    )

    assert len(catalogue) == 6
    assert len({item.provider for item in catalogue}) == 6
    assert len({item.morphology.feature_class for item in catalogue}) == 6
    for item in catalogue:
        density = item.morphology.column_density_m2
        assert density.shape == y_grid.shape
        assert np.all(np.isfinite(density))
        assert np.all(density >= 0.0)
        assert np.isclose(np.max(density), PEAK_DENSITY_M2, rtol=1e-12)
        assert not density.flags.writeable


def test_rotation_and_asymmetry_are_independent_controls(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    aligned = rotated_asymmetric_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        angle_deg=0.0,
        asymmetry=0.0,
    )
    rotated = rotated_asymmetric_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        angle_deg=30.0,
        asymmetry=0.0,
    )
    skewed = rotated_asymmetric_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        angle_deg=0.0,
        asymmetry=0.6,
    )

    aligned_moments = _density_moments(
        aligned.morphology.column_density_m2, y_grid, z_grid
    )
    rotated_moments = _density_moments(
        rotated.morphology.column_density_m2, y_grid, z_grid
    )
    skewed_moments = _density_moments(
        skewed.morphology.column_density_m2, y_grid, z_grid
    )
    assert abs(aligned_moments[4]) < 1e-14
    assert abs(rotated_moments[4]) > 1e-12
    assert abs(aligned_moments[0]) < 1e-12
    assert skewed_moments[0] > 1e-6
    assert rotated.parameters["angle_deg"] == 30.0
    assert skewed.parameters["asymmetry"] == 0.6


def test_transverse_modulation_exposes_period_phase_and_contrast(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    smooth = transverse_modulated_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        modulation_period_um=6.0,
        modulation_contrast=0.0,
    )
    modulated = transverse_modulated_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        modulation_period_um=6.0,
        modulation_contrast=0.7,
        modulation_phase_rad=0.4,
    )

    assert not np.allclose(
        smooth.morphology.column_density_m2,
        modulated.morphology.column_density_m2,
        rtol=1e-6,
        atol=0.0,
    )
    assert modulated.parameters["modulation_period_um"] == 6.0
    assert modulated.parameters["modulation_contrast"] == 0.7
    assert modulated.parameters["modulation_phase_rad"] == 0.4


def test_core_halo_exposes_two_scales_fraction_and_displacement(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    core_only = core_halo_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        core_sigma_y_um=7.0,
        core_sigma_z_um=2.0,
        halo_sigma_y_um=20.0,
        halo_sigma_z_um=7.0,
        halo_integrated_fraction=0.0,
    )
    core_halo = core_halo_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        core_sigma_y_um=7.0,
        core_sigma_z_um=2.0,
        halo_sigma_y_um=20.0,
        halo_sigma_z_um=7.0,
        halo_integrated_fraction=0.30,
        halo_offset_y_um=6.0,
        halo_offset_z_um=2.0,
    )

    core_moments = _density_moments(
        core_only.morphology.column_density_m2, y_grid, z_grid
    )
    mixed_moments = _density_moments(
        core_halo.morphology.column_density_m2, y_grid, z_grid
    )
    assert mixed_moments[0] > core_moments[0]
    assert mixed_moments[1] > core_moments[1]
    assert mixed_moments[2] > 2.0 * core_moments[2]
    assert mixed_moments[3] > 2.0 * core_moments[3]
    assert core_halo.parameters["halo_integrated_fraction"] == 0.30
    assert core_halo.parameters["halo_offset_y_um"] == 6.0
    assert core_halo.morphology.feature_class == "two_scale_core_halo"
    assert "thermodynamic interpretation" in core_halo.morphology.description
    assert not core_halo.morphology.column_density_m2.flags.writeable


@pytest.mark.parametrize("peak_count", [1, 2, 3, 4, 5])
def test_peak_chain_recovers_declared_number_of_resolved_peaks(
    truth_grid: tuple[np.ndarray, np.ndarray],
    peak_count: int,
) -> None:
    y_grid, z_grid = truth_grid
    truth = gaussian_peak_chain(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        peak_count=peak_count,
        separation_um=9.0,
        longitudinal_peak_sigma_um=1.2,
        transverse_peak_sigma_um=1.5,
        peak_contrast=1.0,
    )
    centre_row = int(np.argmin(np.abs(z_grid[:, 0])))
    assert (
        _local_maximum_count(truth.morphology.column_density_m2[centre_row])
        == peak_count
    )
    assert truth.parameters["peak_count"] == peak_count
    assert truth.parameters["separation_um"] == 9.0


def test_peak_width_separation_and_contrast_change_only_declared_inputs(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    baseline = gaussian_peak_chain(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        peak_count=3,
        separation_um=8.0,
        longitudinal_peak_sigma_um=1.5,
        peak_contrast=0.9,
    )
    wider = gaussian_peak_chain(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        peak_count=3,
        separation_um=8.0,
        longitudinal_peak_sigma_um=3.0,
        peak_contrast=0.9,
    )
    separated = gaussian_peak_chain(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        peak_count=3,
        separation_um=11.0,
        longitudinal_peak_sigma_um=1.5,
        peak_contrast=0.9,
    )
    lower_contrast = gaussian_peak_chain(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        peak_count=3,
        separation_um=8.0,
        longitudinal_peak_sigma_um=1.5,
        peak_contrast=0.3,
    )

    baseline_density = baseline.morphology.column_density_m2
    assert not np.array_equal(baseline_density, wider.morphology.column_density_m2)
    assert not np.array_equal(baseline_density, separated.morphology.column_density_m2)
    assert not np.array_equal(
        baseline_density, lower_contrast.morphology.column_density_m2
    )
    assert wider.parameters["separation_um"] == baseline.parameters["separation_um"]
    assert separated.parameters["longitudinal_peak_sigma_um"] == baseline.parameters[
        "longitudinal_peak_sigma_um"
    ]
    assert lower_contrast.parameters["peak_contrast"] == 0.3


def test_hole_depth_and_notch_depth_are_explicit_controls(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    no_hole = central_hole_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        hole_contrast=0.0,
    )
    hole = central_hole_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        hole_contrast=1.0,
    )
    centre = np.unravel_index(np.argmin(y_grid**2 + z_grid**2), y_grid.shape)
    assert no_hole.morphology.column_density_m2[centre] > 0.99 * PEAK_DENSITY_M2
    assert hole.morphology.column_density_m2[centre] == pytest.approx(0.0, abs=1e-6)

    smooth_curve = curved_notched_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        notch_centre_y_um=3.0,
        notch_contrast=0.0,
    )
    notched_curve = curved_notched_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        notch_centre_y_um=3.0,
        notch_contrast=0.8,
    )
    target = np.unravel_index(
        np.argmin((y_grid * 1e6 - 3.0) ** 2 + (z_grid * 1e6 - 0.0675) ** 2),
        y_grid.shape,
    )
    assert (
        notched_curve.morphology.column_density_m2[target]
        < 0.3 * smooth_curve.morphology.column_density_m2[target]
    )


def test_expanded_shifted_input_reaches_field_boundary(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    truth = expanded_shifted_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        centre_y_um=28.0,
        centre_z_um=4.0,
        longitudinal_sigma_um=25.0,
        transverse_sigma_um=9.0,
    )
    density = truth.morphology.column_density_m2
    centre_y, centre_z, *_ = _density_moments(density, y_grid, z_grid)
    assert centre_y > 10e-6
    assert centre_z > 1e-6
    assert np.max(density[:, -1]) > 0.1 * np.max(density)
    assert truth.morphology.feature_class == "support_challenge"


def test_latent_phase_pair_has_identical_density_and_distinct_metadata(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    base = central_hole_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
    )
    first, second = latent_phase_impossibility_pair(
        base,
        phase_labels=("zero_winding", "unit_winding"),
        equivalence_group="same_density_different_phase",
    )

    assert np.array_equal(
        first.morphology.column_density_m2,
        second.morphology.column_density_m2,
    )
    assert first.morphology.column_density_m2 is second.morphology.column_density_m2
    assert first.latent_metadata["latent_phase_label"] == "zero_winding"
    assert second.latent_metadata["latent_phase_label"] == "unit_winding"
    assert (
        first.observational_equivalence_group
        == second.observational_equivalence_group
        == "same_density_different_phase"
    )


def test_controlled_sequence_applies_declared_schedules_without_dynamics(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    base = rotated_asymmetric_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        longitudinal_sigma_um=8.0,
        transverse_sigma_um=2.5,
        angle_deg=0.0,
        asymmetry=0.0,
    )
    sequence = controlled_unknown_sequence(
        y_grid,
        z_grid,
        base=base,
        times_s=(0.0, 0.01, 0.02),
        amplitude_scales=(1.0, 0.8, 0.6),
        translation_y_um=(0.0, 2.0, 4.0),
    )

    assert isinstance(sequence, AnalyticTruthSequence)
    assert sequence.times_s == (0.0, 0.01, 0.02)
    assert sequence.parameters["dynamics_model"] == "none"
    assert np.allclose(
        sequence.frames[0].morphology.column_density_m2,
        base.morphology.column_density_m2,
        rtol=1e-14,
        atol=0.0,
    )
    totals = [
        float(np.sum(frame.morphology.column_density_m2))
        for frame in sequence.frames
    ]
    centres = [
        _density_moments(
            frame.morphology.column_density_m2,
            y_grid,
            z_grid,
        )[0]
        for frame in sequence.frames
    ]
    assert totals[1] / totals[0] == pytest.approx(0.8, rel=2e-3)
    assert totals[2] / totals[0] == pytest.approx(0.6, rel=2e-3)
    assert centres[1] - centres[0] == pytest.approx(2.0e-6, abs=2e-8)
    assert centres[2] - centres[0] == pytest.approx(4.0e-6, abs=2e-8)
    assert all(
        not frame.morphology.column_density_m2.flags.writeable
        for frame in sequence.frames
    )


def test_controlled_sequence_scale_and_abrupt_peak_are_independent(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    base = rotated_asymmetric_cloud(
        y_grid,
        z_grid,
        peak_column_density_m2=PEAK_DENSITY_M2,
        longitudinal_sigma_um=6.0,
        transverse_sigma_um=2.0,
        angle_deg=0.0,
        asymmetry=0.0,
    )
    sequence = controlled_unknown_sequence(
        y_grid,
        z_grid,
        base=base,
        times_s=(0.0, 0.01, 0.02),
        scale_y=(1.0, 1.4, 1.4),
        sudden_peak_frame_index=2,
        sudden_peak_relative_amplitude=0.5,
        sudden_peak_centre_y_um=20.0,
        sudden_peak_centre_z_um=0.0,
        sudden_peak_sigma_y_um=1.5,
        sudden_peak_sigma_z_um=1.5,
    )

    first_variance = _density_moments(
        sequence.frames[0].morphology.column_density_m2, y_grid, z_grid
    )[2]
    second_variance = _density_moments(
        sequence.frames[1].morphology.column_density_m2, y_grid, z_grid
    )[2]
    target = np.unravel_index(
        np.argmin((y_grid * 1e6 - 20.0) ** 2 + (z_grid * 1e6) ** 2),
        y_grid.shape,
    )
    before = sequence.frames[1].morphology.column_density_m2[target]
    after = sequence.frames[2].morphology.column_density_m2[target]
    assert second_variance > 1.8 * first_variance
    assert after - before > 0.45 * PEAK_DENSITY_M2
    assert sequence.frames[1].parameters["sudden_peak_active"] is False
    assert sequence.frames[2].parameters["sudden_peak_active"] is True


def test_controlled_sequence_copies_raw_density_input(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    raw_density = np.exp(
        -0.5 * ((y_grid / 8e-6) ** 2 + (z_grid / 3e-6) ** 2)
    )
    sequence = controlled_unknown_sequence(
        y_grid,
        z_grid,
        base=raw_density,
        times_s=(0.0, 0.01),
    )
    initial_copy = np.array(sequence.frames[0].morphology.column_density_m2)
    raw_density[:] = 0.0
    assert np.array_equal(
        sequence.frames[0].morphology.column_density_m2, initial_copy
    )
    assert sequence.parameters["source_provider"] == "array_input"


def test_provider_parameter_validation(
    truth_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    y_grid, z_grid = truth_grid
    with pytest.raises(ValueError, match="between one and five"):
        gaussian_peak_chain(
            y_grid,
            z_grid,
            peak_column_density_m2=PEAK_DENSITY_M2,
            peak_count=6,
        )
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        transverse_modulated_cloud(
            y_grid,
            z_grid,
            peak_column_density_m2=PEAK_DENSITY_M2,
            modulation_contrast=1.1,
        )
    with pytest.raises(ValueError, match="match peak_count"):
        gaussian_peak_chain(
            y_grid,
            z_grid,
            peak_column_density_m2=PEAK_DENSITY_M2,
            peak_count=3,
            relative_peak_weights=(1.0, 0.8),
        )
    with pytest.raises(ValueError, match="strictly increasing"):
        controlled_unknown_sequence(
            y_grid,
            z_grid,
            base=np.ones_like(y_grid),
            times_s=(0.0, 0.0),
        )
    with pytest.raises(ValueError, match="one value per sequence frame"):
        controlled_unknown_sequence(
            y_grid,
            z_grid,
            base=np.ones_like(y_grid),
            times_s=(0.0, 0.1),
            amplitude_scales=(1.0,),
        )
    with pytest.raises(ValueError, match="select a sequence frame"):
        controlled_unknown_sequence(
            y_grid,
            z_grid,
            base=np.ones_like(y_grid),
            times_s=(0.0, 0.1),
            sudden_peak_frame_index=2,
        )
