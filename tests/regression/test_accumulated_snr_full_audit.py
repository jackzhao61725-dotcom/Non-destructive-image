import json
from pathlib import Path

import pytest

from scripts.audit_accumulated_snr_physics_and_code import run_audit


ROOT = Path(__file__).resolve().parents[2]


def test_full_accumulated_snr_audit_outputs_and_verdict(tmp_path: Path) -> None:
    outputs = run_audit(ROOT / "configs" / "dissertation_plots_v1.json")
    assert all(path.exists() for path in outputs.values())
    with outputs["summary"].open(encoding="utf-8") as handle:
        summary = json.load(handle)

    assert summary["verdict"].startswith("B.")
    exponents = {row["range"]: row for row in summary["scaling_exponents"]}
    assert exponents["full_0.75_to_5.0_GHz"]["clean_loss_nmax_exponent"] == pytest.approx(2.0, abs=0.01)
    assert exponents["full_0.75_to_5.0_GHz"]["heating_reabs_nmax_exponent"] > 1.9
    assert len(summary["critical_or_quantitative_issues"]) == 2
