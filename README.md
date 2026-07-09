# Non-destructive Imaging Simulator for Ultracold 166Er Condensates

This repository contains a notebook-equivalent simulator and
calibration-aware optimisation framework for continuous non-destructive imaging
of ultracold `166Er` condensates.

The original notebook remains the scientific reference implementation:

```text
1 calculations revised 2  multishot  6  extended.ipynb
```

The helper package in `src/non_destructive_image/` preserves
notebook-equivalent behaviour while making the simulator easier to test,
document, and extend.

## Current Version 1 Status

Version 1 is closed for the migrated simulator core:

- Atomic Model helpers are implemented and tested.
- Light-Atom Interaction helpers are implemented and tested.
- Imaging helpers are implemented and tested, including PCI, DGI, and Faraday.
- Camera helpers are implemented and tested, including deterministic and
  stochastic camera paths.
- Deterministic multi-shot sequence helpers are implemented and tested.
- Deterministic Faraday optimisation helpers are implemented and tested.
- Absorption / RAI calibration-readiness helpers are implemented and tested.
- Representative dissertation outputs have been generated under
  `results/faraday_optimisation_v1/`.

## Architecture Overview

```text
Atomic Model
  -> Light-Atom Interaction
  -> Imaging
  -> Camera
  -> Multi-shot
  -> Analysis / Calibration
```

The optimisation and calibration layers sit above the notebook-equivalent
simulator core. They consume validated helper outputs rather than changing the
underlying physics.

## Key Repository Locations

- `src/non_destructive_image/` - simulator helper package.
- `tests/` - regression and helper tests.
- `notebook_sections/` - exported notebook sections for audit and validation.
- `scripts/` - validation, baseline, bundle, and result-generation scripts.
- `configs/` - configurable result-generation inputs.
- `results/` - representative dissertation output tables and figures.
- `regression/` - regression baselines used by tests.
- `docs/` - architecture, migration, optimisation, calibration, and
  dissertation-facing documentation.

## Validation

On Windows PowerShell:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

## Regenerate Representative Results

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
python scripts\generate_dissertation_results.py --config configs\dissertation_results_v1.json
```

The generated outputs are written to:

```text
results/faraday_optimisation_v1/
```

## Current Limitations

- The generated outputs are Version 1 representative / uncalibrated results.
- `kappa_F = 1.0` remains a placeholder convention.
- No experimental RAI / absorption calibration has been applied yet.
- No `kappa_F` fitting has been implemented yet.
- No microscopic Faraday model has been introduced.
- No full multi-parameter optimisation has been implemented yet.

## Branches

- `dissertation-v1-snapshot` is the stable dissertation snapshot branch.
- `clean-v1-working` is the cleaned working branch intended for continuing MSc
  project work.
