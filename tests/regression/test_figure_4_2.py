from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.generate_figure_4_2 import DEFAULT_CONFIG, generate


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, Path], dict]:
    output = tmp_path_factory.mktemp("figure_4_2")
    outputs = generate(DEFAULT_CONFIG, output)
    values = json.loads(outputs["values"].read_text(encoding="utf-8"))
    return outputs, values


def test_figure_4_2_canonical_gate_passes(generated: tuple[dict[str, Path], dict]) -> None:
    _, values = generated
    assert values["canonical_gate"]["passed"] is True
    assert values["optical_input"]["phase_peak_rad"] == pytest.approx(0.20294165287929014)
    assert values["optical_input"]["kappa_F"] == pytest.approx(-45 / 91)
    assert values["optical_input"]["numerical_aperture"] == pytest.approx(0.13)
    assert values["optical_input"]["theta_f_at_scalar_phase_peak_rad"] == pytest.approx(
        -0.1003557624128358
    )
    assert values["optical_input"]["scattered_power_throughput"] == pytest.approx(
        0.9302787994283651
    )
    assert values["readouts"]["pci"]["stats"]["max"] == pytest.approx(1.2780199161462586)
    assert values["readouts"]["dgi"]["stats"]["max"] == pytest.approx(0.03257013971289065)
    assert values["readouts"]["dark_field_faraday"]["stats"]["max"] == pytest.approx(
        0.008018767598998892
    )
    dual_port = values["readouts"]["dual_port_faraday"]
    assert dual_port["signal_stats"]["centre_value"] == pytest.approx(-0.17833439662712702)
    assert dual_port["signal_stats"]["min"] == pytest.approx(-0.17833439662712702)
    assert dual_port["port_H_stats"]["centre_value"] == pytest.approx(0.41102438482000514)
    assert dual_port["port_V_stats"]["centre_value"] == pytest.approx(0.5894419438976197)
    assert dual_port["signal_definition"] == "S = (I_H - I_V)/(I_H + I_V)"
    assert dual_port["port_mapping"] == (
        "I_H = historical notebook I_v; I_V = historical notebook I_u"
    )
    centre_field = values["optical_input"]["centre_filtered_scattered_field"]
    assert centre_field["real"] == pytest.approx(-0.015439786478921252)
    assert centre_field["imag"] == pytest.approx(0.1803899898440997)
    assert centre_field["magnitude_squared"] == pytest.approx(0.03277893544246907)


def test_figure_4_2_outputs_are_thesis_ready(generated: tuple[dict[str, Path], dict]) -> None:
    outputs, _ = generated
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())
    svg = outputs["svg"].read_text(encoding="utf-8")
    assert "(a) PCI" in svg
    assert "(b) DGI" in svg
    assert "(c) Dark-field Faraday" in svg
    assert "(d) Dual-port Faraday" in svg
    assert "Stage" not in svg

    metadata = json.loads(outputs["metadata"].read_text(encoding="utf-8"))
    assert metadata["camera_noise_included"] is False
    assert metadata["multishot_evolution_included"] is False
    assert metadata["port_mapping"] == "I_H = historical notebook I_v; I_V = historical notebook I_u"
    assert "kappa_F=-45/91" in metadata["atomic_response_status"]
    assert "remain to be measured" in metadata["apparatus_calibration_status"]
