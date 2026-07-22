from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.generate_figure_5_1 import DEFAULT_CONFIG, generate


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, Path], dict]:
    output = tmp_path_factory.mktemp("figure_5_1")
    outputs = generate(DEFAULT_CONFIG, output)
    values = json.loads(outputs["values"].read_text(encoding="utf-8"))
    return outputs, values


def test_figure_5_1_canonical_context(generated: tuple[dict[str, Path], dict]) -> None:
    _, values = generated
    reference = values["canonical_gate"]["reference"]
    assert values["canonical_gate"]["passed"] is True
    assert reference["detuning_ghz"] == pytest.approx(1.5)
    assert reference["fluence_mw_us"] == pytest.approx(300.0)
    assert reference["phase_peak_rad"] == pytest.approx(0.20294165287929014)
    assert reference["photoelectrons_per_incident_I0_pixel"] == pytest.approx(735.2695759176155)
    assert reference["snr"] == pytest.approx(
        {
            "Faraday dark-field": 7.800030869978132,
            "Faraday dual-port": 15.34463436369928,
        }
    )
    assert values["fixed_context"]["central_block_shape"] == [5, 5]
    assert values["fixed_context"]["camera_output_shape"] == [153, 153]
    assert values["fixed_context"]["effective_magnification"] == pytest.approx(10.0)
    assert values["fixed_context"]["kappa_F"] == pytest.approx(-45 / 91)
    assert values["fixed_context"]["numerical_aperture"] == pytest.approx(0.13)
    assert values["fixed_context"]["read_noise_electrons_per_pixel_per_port"] == pytest.approx(
        0.7
    )
    assert values["selected_frames"]["fluence_mw_us"] == pytest.approx([90.0, 150.0, 300.0])
    assert values["selected_frames"]["snr"] == pytest.approx(
        {
            "Faraday dark-field": [
                3.6877291730915167,
                5.150743161494631,
                7.800030869978132,
            ],
            "Faraday dual-port": [
                8.391572042384801,
                10.843075956189569,
                15.34463436369928,
            ],
        }
    )


def test_figure_5_1_outputs_are_complete(generated: tuple[dict[str, Path], dict]) -> None:
    outputs, values = generated
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())
    svg = outputs["svg"].read_text(encoding="utf-8")
    assert "Fixed detuning" in svg
    assert "Dark-field" in svg
    assert "Dual-port" in svg
    assert "(a) PCI" not in svg
    metadata = json.loads(outputs["metadata"].read_text(encoding="utf-8"))
    assert "-45/91" in metadata["faraday_boundary"]
    assert "apparatus-level response remains to be calibrated" in metadata["faraday_boundary"]
    for mode, limits in values["scan"]["snr_ranges"].items():
        assert mode in {"Faraday dark-field", "Faraday dual-port"}
        assert 0 < limits["min"] < limits["max"]
