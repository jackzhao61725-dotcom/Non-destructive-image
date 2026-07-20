import csv
import json
from pathlib import Path


from scripts.run_performance_validation import (
    REQUIRED_LEDGER_COLUMNS,
    evaluate_canonical_gate,
)


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "configs" / "performance_validation_v1.json"
OUTPUT_DIR = ROOT / "results" / "performance_validation_v1"


def test_canonical_performance_gate_matches_maintained_results() -> None:
    with CONFIG_PATH.open(encoding="utf-8") as handle:
        config = json.load(handle)
    rows, checks, invariants = evaluate_canonical_gate(
        config,
        CONFIG_PATH,
        OUTPUT_DIR / "canonical_gate.csv",
    )

    assert len(rows) == 8
    assert all(check["passed"] for check in checks)
    assert invariants == {
        "accepted_frames": 10,
        "post_sequence_loss": 0.2806023117605313,
        "next_pulse_loss": 0.30829044050175103,
        "threshold_crossing_state_included": "False",
    }
    faraday_rows = [row for row in rows if row["mode"].startswith("Faraday")]
    assert len(faraday_rows) == 4
    assert all("uncalibrated" in row["calibration status"] for row in faraday_rows)
    assert all(row["kappa_F"] == "1.0 placeholder" for row in faraday_rows)


def test_canonical_gate_evidence_has_complete_provenance() -> None:
    with (OUTPUT_DIR / "canonical_gate.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    with (OUTPUT_DIR / "canonical_gate.json").open(encoding="utf-8") as handle:
        evidence = json.load(handle)

    assert rows
    assert set(REQUIRED_LEDGER_COLUMNS).issubset(rows[0])
    assert all(all(row[column] for column in REQUIRED_LEDGER_COLUMNS) for row in rows)
    assert evidence["passed"] is True
    assert evidence["faraday_boundary"]["calibrated_absolute_prediction"] is False
    assert evidence["source_sha256"]
