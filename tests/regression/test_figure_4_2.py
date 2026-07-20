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
    assert values["optical_input"]["scattered_power_throughput"] == pytest.approx(0.7115304533861209)
    assert values["readouts"]["pci"]["stats"]["max"] == pytest.approx(1.1585677695076864)
    assert values["readouts"]["dgi"]["stats"]["max"] == pytest.approx(0.01595651906354294)
    assert values["readouts"]["dark_field_faraday"]["stats"]["max"] == pytest.approx(
        0.015956443479481053
    )
    assert values["readouts"]["dual_port_faraday"]["signal_stats"]["max"] == pytest.approx(
        0.25116900542198173
    )
    centre_field = values["optical_input"]["centre_filtered_scattered_field"]
    assert centre_field["real"] == pytest.approx(-0.010274925557004756)
    assert centre_field["imag"] == pytest.approx(0.12631881680684415)
    assert centre_field["magnitude_squared"] == pytest.approx(0.016062017574683046)


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
    assert "uncalibrated" in metadata["calibration_status"]
