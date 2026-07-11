import json
from pathlib import Path

import numpy as np

from scripts.generate_full_multishot_accumulated_snr import matched_filter_snr


ROOT = Path(__file__).resolve().parents[2]


def test_matched_filter_snr_orders_shot_and_read_noise() -> None:
    image = np.array([[1.0, 1.1], [0.9, 1.2]])
    roi = np.ones_like(image, dtype=bool)
    result = matched_filter_snr(image, 1.0, 1000.0, 7.0, roi)
    assert result["snr_shot_plus_read_noise"] <= result["snr_shot_noise_only"]
    assert result["signal_l2_e"] > 0
    assert np.isfinite(result["shot_plus_read_effective_e"])


def test_generated_full_multishot_summary_contract() -> None:
    output = ROOT / "results" / "dissertation_plots_v1" / "full_multishot_accumulated_snr"
    required = [
        "full_multishot_accumulated_snr.svg",
        "full_multishot_accumulated_snr_data.csv",
        "framewise_snr_sequences.csv",
        "model_comparison.csv",
        "full_multishot_accumulated_snr_summary.json",
        "metadata.json",
    ]
    assert all((output / name).exists() for name in required)
    with (output / "full_multishot_accumulated_snr_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)
    assert summary["checks"]["passed"] is True
    assert summary["checks"]["pci_shot_plus_read_le_shot_all"] is True
    assert summary["checks"]["dgi_shot_plus_read_le_shot_all"] is True
    assert summary["verdict"].startswith("A.")
