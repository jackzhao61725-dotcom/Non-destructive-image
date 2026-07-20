# Performance Validation V1

## Current stage

Canonical gate: **PASS**.
Branch: `work/full-performance-validation`. Source commit: `bbcc8d93293634ecee7fc4c2ea13ef5571b129ff`.

This stage verifies the maintained fixed 228-pixel matched-ROI, evolving-cloud, heating + initial-density reabsorption sequence. Peak-pixel, resolution-element, and other estimators are not mixed into this table.
The camera model is `Thorlabs DCC3260M`. QE=0.6 and read noise 3 e- rms per pixel and readout are provisional dissertation screening values; neither is a measurement of the installed camera.

## Validation commands

- Prerequisites intentionally skipped; numerical contract only was evaluated.

## Canonical four-mode results

| Mode | Noise | Actual | Expected | Absolute difference |
| --- | --- | ---: | ---: | ---: |
| PCI | shot_noise_only | 104.098235973135 | 104.098235973135 | 0 |
| PCI | shot_plus_read_noise | 103.668651367624 | 103.668651367624 | 0 |
| DGI | shot_noise_only | 54.3619567765478 | 54.3619567765478 | 0 |
| DGI | shot_plus_read_noise | 34.3032571045613 | 34.3032571045613 | 0 |
| Faraday dark-field | shot_noise_only | 55.7783373129339 | 55.7783373129339 | 0 |
| Faraday dark-field | shot_plus_read_noise | 34.708770831391 | 34.708770831391 | 0 |
| Faraday dual-port | shot_noise_only | 107.481977345106 | 107.481977345106 | 0 |
| Faraday dual-port | shot_plus_read_noise | 106.555526352218 | 106.555526352218 | 0 |

Accepted frames: `10` (indices `0..9`). The threshold-crossing pulse is excluded.
Post-sequence loss: `0.2806023117605313`.
Next-pulse loss: `0.30829044050175103`.

## Calibration boundary

PCI and DGI rows are Version 1 notebook-aligned numerical results. Faraday rows use the phenomenological placeholder `kappa_F=1.0`; they are uncalibrated structural comparisons, not calibrated absolute or experimental predictions.

## Issue register

- None at the canonical gate.

## Evidence

Machine-readable evidence is in `results/performance_validation_v1/canonical_gate.csv` and `canonical_gate.json`. The source full-multishot CSV, Faraday ledger, metadata, hashes, command output, and exact tolerances are recorded there.

Later scientific, statistical, numerical, computational, and software validation stages remain separate; this report does not collapse them into a single performance score.
