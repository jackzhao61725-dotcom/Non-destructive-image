from __future__ import annotations

import csv
import json
from pathlib import Path


OUTPUT_DIR = Path("results/linear_approximation_audit")


def test_linear_approximation_audit_outputs_exist() -> None:
    expected = [
        OUTPUT_DIR / "linear_approximation_summary.json",
        OUTPUT_DIR / "phase_rotation_ranges.csv",
        OUTPUT_DIR / "small_angle_error_table.csv",
        OUTPUT_DIR / "detuning_scaling_error_table.csv",
        OUTPUT_DIR / "metadata.json",
        Path("docs/linear_approximation_validity_audit.md"),
    ]
    for path in expected:
        assert path.exists(), path


def test_linear_approximation_audit_records_finite_phase_and_rotation_ranges() -> None:
    summary = json.loads((OUTPUT_DIR / "linear_approximation_summary.json").read_text(encoding="utf-8"))
    ranges = summary["canonical_ranges"]
    assert ranges["max_abs_phi_rad"] > 0
    assert ranges["central_abs_phi_rad"] > 0
    assert ranges["max_abs_theta_rad"] > 0
    assert ranges["central_abs_theta_rad"] > 0
    assert ranges["max_abs_phi_rad"] == ranges["max_abs_theta_rad"]


def test_linear_approximation_audit_tables_have_required_columns() -> None:
    with (OUTPUT_DIR / "small_angle_error_table.csv").open(newline="", encoding="utf-8") as handle:
        small_angle_rows = list(csv.DictReader(handle))
    assert small_angle_rows
    required_small_angle_columns = {
        "quantity",
        "sample",
        "value_rad",
        "exp_i_phi_vs_1_plus_i_phi_relative_field_error",
        "dual_port_sin_2theta_vs_2theta_relative_error",
    }
    assert required_small_angle_columns.issubset(set(small_angle_rows[0]))

    with (OUTPUT_DIR / "detuning_scaling_error_table.csv").open(newline="", encoding="utf-8") as handle:
        detuning_rows = list(csv.DictReader(handle))
    assert detuning_rows
    required_detuning_columns = {
        "detuning_hz",
        "dimensionless_delta",
        "phase_scaling_relative_error",
        "scattering_scaling_relative_error",
    }
    assert required_detuning_columns.issubset(set(detuning_rows[0]))


def test_linear_approximation_audit_metadata_records_provenance() -> None:
    metadata = json.loads((OUTPUT_DIR / "metadata.json").read_text(encoding="utf-8"))
    assert "configs\\notebook_v1_defaults.json" in metadata["config_files_used"] or (
        "configs/notebook_v1_defaults.json" in metadata["config_files_used"]
    )
    assert metadata["simulator_physics_changed"] is False
    assert metadata["helper_apis_changed"] is False
    assert metadata["notebook_logic_changed"] is False
