# Reproducibility

## Scope

This repository is prepared as a reproducible Version 1 simulation and
dissertation-output workflow. The goal is to make the current notebook-aligned
figures, representative plots, and validation checks rerunnable before a future
archived Zenodo release.

This does not create a GitHub release, does not trigger Zenodo deposition, and
does not claim an experimental calibration.

## Environment

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
```

## Validation

Run:

```powershell
pytest -q
python scripts\validate_notebook_sections.py
```

The notebook-section validator checks that unmigrated exported notebook
sections remain in sync and that migrated sections are syntactically valid.

## One-Command Output Reproduction

Run:

```powershell
python scripts\run_all_dissertation_figures.py
```

The script runs approved output-generation scripts in a controlled order and
writes:

```text
results/reproducibility_manifest.json
```

The manifest records:

- UTC timestamp;
- git commit hash;
- Python executable and version;
- scripts run;
- configs used;
- expected output files;
- pending figure workflows not yet merged into `main`.

## Included Generators

The current `run_all` entry point includes:

- notebook-section validation;
- condensate recovery;
- scalar phase recovery;
- PCI recovery;
- DGI recovery;
- Faraday recovery;
- deterministic camera recovery;
- stochastic camera recovery;
- deterministic multishot recovery;
- noisy PCI multishot filmstrip recovery;
- condensate three-view projection;
- representative Faraday optimisation result generation;
- detuning trade-off physics plot generation;
- finite-phase / finite-Faraday-rotation linear-approximation audit.

## Parameter Provenance

No figure without parameter provenance.

Every figure or plot should record:

- config files used;
- physical parameters;
- optical parameters;
- camera/noise parameters where relevant;
- sweep parameters where relevant;
- units;
- normalisation rules;
- calibration status.

Notebook-aligned recovery outputs use:

```text
configs/notebook_v1_defaults.json
```

Representative Faraday optimisation outputs use:

```text
configs/dissertation_results_v1.json
```

The dissertation detuning trade-off plot uses:

```text
configs/dissertation_plots_v1.json
```

The linear-approximation audit uses the notebook defaults and records whether
the dissertation plot config is present:

```text
configs/notebook_v1_defaults.json
configs/dissertation_plots_v1.json
```

## Current Scientific Boundary

The current outputs are Version 1 representative / uncalibrated results.

Important limitations:

- `kappa_F = 1.0` remains a phenomenological placeholder;
- no experimental RAI / absorption calibration has been applied;
- no `kappa_F` fitting has been implemented;
- no final calibrated operating-point prediction is claimed;
- future calibrated outputs should use explicit calibration config files.

## Before Zenodo Archival

Before creating a GitHub release for Zenodo:

1. run full validation;
2. run `scripts\run_all_dissertation_figures.py`;
3. inspect `results/reproducibility_manifest.json`;
4. verify README, `CITATION.cff`, `.zenodo.json`, and license status;
5. choose and add a license if public reuse is intended.
