"""Pure tabular summaries shared by benchmark generation and reporting."""

from __future__ import annotations

from typing import Any

import numpy as np


def _successful(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def aggregate_held_out_trials(
    rows: tuple[dict[str, Any], ...] | list[dict[str, Any]],
) -> tuple[dict[str, Any], ...]:
    """Summarise held-out trials while retaining failures in the denominator."""

    output: list[dict[str, Any]] = []
    for readout in ("dual_port", "dark_field"):
        fluences = sorted(
            {
                float(row["fluence_mw_us"])
                for row in rows
                if row["readout"] == readout
            }
        )
        for fluence in fluences:
            trials = [
                row
                for row in rows
                if row["readout"] == readout
                and float(row["fluence_mw_us"]) == fluence
            ]
            successful = [row for row in trials if _successful(row["success"])]
            common = {
                "readout": readout,
                "fluence_mw_us": fluence,
                "trial_count": len(trials),
                "successful_trials": len(successful),
                "success_fraction": len(successful) / len(trials),
            }
            if not successful:
                output.append(
                    {
                        **common,
                        "median_supported_band_relative_l2_error": "",
                        "minimum_supported_band_relative_l2_error": "",
                        "maximum_supported_band_relative_l2_error": "",
                        "median_absolute_integrated_density_relative_error": "",
                    }
                )
                continue
            supported = np.asarray(
                [float(row["supported_band_relative_l2_error"]) for row in successful]
            )
            integrated = np.abs(
                np.asarray(
                    [float(row["integrated_density_relative_error"]) for row in successful]
                )
            )
            output.append(
                {
                    **common,
                    "median_supported_band_relative_l2_error": float(
                        np.median(supported)
                    ),
                    "minimum_supported_band_relative_l2_error": float(
                        np.min(supported)
                    ),
                    "maximum_supported_band_relative_l2_error": float(
                        np.max(supported)
                    ),
                    "median_absolute_integrated_density_relative_error": float(
                        np.median(integrated)
                    ),
                }
            )
    return tuple(output)
