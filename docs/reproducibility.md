# Reproducibility

## Environment

```powershell
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
```

## Validation

Run the full test suite:

```powershell
pytest -q
```

Run the canonical gate against the maintained generated data:

```powershell
python scripts\run_performance_validation.py --skip-prerequisites
```

The gate verifies the four accumulated matched-template SNR pairs, strict
10-frame stopping, the 228-pixel common ROI and the uncalibrated Faraday
boundary.

## One-command regeneration

Inspect the ordered plan:

```powershell
python scripts\run_all_dissertation_figures.py --dry-run
```

Regenerate all approved historical recovery outputs, current dissertation
figures, numerical audits and the canonical gate:

```powershell
python scripts\run_all_dissertation_figures.py
```

The script writes `results/reproducibility_manifest.json`, containing the
Python environment, Git commit, commands, configs and expected outputs.

The complete run includes the expensive Figure 5.2 detuning-fluence screen.
Allow several minutes.

## Active configs

| Role | Config |
| --- | --- |
| historical notebook regression | `configs/notebook_v1_defaults.json` |
| active optical/detector/sequence reference | `configs/dissertation_v2_dcc3260m.json` |
| shared screening and accumulated-SNR scans | `configs/dissertation_plots_v1.json` |
| figures | `configs/figure_4_2.json`, `figure_5_1.json`, `figure_5_2.json`, `figure_5_4.json` |
| thesis numerical audit | `configs/thesis_numerical_contract_v1.json` |
| canonical gate | `configs/performance_validation_v1.json` |

## Dependency order

The accumulated-SNR scaling table is generated before the full multishot
comparison because the latter imports the former as an explicitly labelled
clean-loss reference. The active CSV column is `SNR_acc`.

Figure 5.4 uses the same physical and image-quality contract as Figure 5.2.
Changing magnification, QE, read noise, numerical aperture, `kappa_F` or the
disturbance model therefore requires regeneration of the dependent Chapter 5
outputs and the canonical gate.

## Focused commands

```powershell
python scripts\generate_accumulated_snr_invariance_plot.py --config configs\dissertation_plots_v1.json
python scripts\generate_full_multishot_accumulated_snr.py --config configs\dissertation_plots_v1.json
python scripts\generate_figure_4_2.py --config configs\figure_4_2.json
python scripts\generate_figure_5_1.py --config configs\figure_5_1.json
python scripts\generate_figure_5_2.py --config configs\figure_5_2.json
python scripts\generate_figure_5_4_snr_panel.py --config configs\figure_5_4.json
python scripts\audit_thesis_numerical_consistency.py --config configs\thesis_numerical_contract_v1.json
```

## Provenance requirements

Every retained output must record:

- source config paths;
- detuning, power, exposure and imaging axis;
- optical sampling and detector parameters;
- signal normalisation and estimator/ROI;
- sequence and stopping model;
- random seed when a stochastic frame is rendered;
- calibration status;
- output path and Git commit.

Generated camera images are stochastic illustrations. Reported SNR is computed
from expected counts and the analytic noise variance.

## Scientific boundary

The current result set is internally reproducible but experimentally
uncalibrated. In particular, `kappa_F=1` does not make the Faraday values
calibrated predictions. Commissioning data should first replace provisional
inputs; validation requires comparison with separate data not used for that
replacement.
