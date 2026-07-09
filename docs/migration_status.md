# Migration Status

This document records the current Version 1 migration status.

The original notebook remains the authoritative scientific implementation. The
helper package is a conservative, regression-tested support layer that preserves
notebook-equivalent behaviour for the migrated core. No physics redesign was
introduced as part of Version 1 migration.

## Version 1 Migrated Core Status

The repository now contains a closed Version 1 migrated simulator core:

- Atomic Model
- Light-Atom Interaction
- Imaging
- Camera
- Stochastic Camera Noise
- Deterministic Multi-shot Core

The repository also contains a deterministic single-variable Faraday
optimisation layer above the migrated core. This layer supports operating-point
evaluation, one-dimensional detuning / probe-power / exposure-time sweeps, and
small summary dictionaries for report use.

The exported notebook sections remain in `notebook_sections/` for auditability,
validation, and future migration work. They are not replaced as the scientific
source of truth.

## Dependency Graph

```text
Original Notebook (authoritative scientific implementation)
  |
  | exported and validated
  v
Notebook Section Exports
  |
  | conservative helper extraction
  v
src/non_destructive_image helper package
  |
  | regression tests and stored baselines
  v
Version 1 migrated core
```

## Helper Validation Status

| Layer | Module / location | Status | Notes |
|---|---|---|---|
| Atomic Model | `src/non_destructive_image/atomic_model.py`, `profiles.py` | Version 1 core migrated and tested | Thomas-Fermi state construction, recoil quantities, and 2D profile helper are available. |
| Light-Atom Interaction | `src/non_destructive_image/light_atom.py` | Version 1 core migrated and tested | Scalar phase, residual OD, intensity, scattering, Faraday scaling, and reabsorption helpers are available. |
| Shared Fourier Core | `src/non_destructive_image/fourier.py`, `imaging.py` | Version 1 core migrated and tested | Shared FFT/pupil propagation and coherent Fourier image recombination are available. |
| PCI orchestration | `src/non_destructive_image/imaging.py` | Version 1 core migrated and tested | `simulate_pci_image(...)` is tested against the PCI/DGI imaging baseline. |
| DGI orchestration | `src/non_destructive_image/imaging.py` | Version 1 core migrated and tested | `simulate_dgi_image(...)` is tested against the PCI/DGI imaging baseline. |
| Faraday orchestration | `src/non_destructive_image/imaging.py` | Version 1 core migrated and tested | `simulate_faraday_image(...)` is tested against the Faraday imaging baseline and preserves the current `kappa_F = 1.0` phenomenological convention. |
| Camera deterministic pipeline | `src/non_destructive_image/camera.py` | Version 1 core migrated and tested | Binning and deterministic normalisation are helperized. |
| Stochastic camera noise | `src/non_destructive_image/camera.py` | Version 1 core migrated and tested | Noise helper uses explicit `np.random.Generator`; no hidden global seed is introduced. |
| Deterministic multi-shot core | `src/non_destructive_image/multishot.py` | Version 1 core migrated and tested | Heating / clean-loss bookkeeping and RMS accumulated SNR are helperized. |
| Faraday optimisation | `src/non_destructive_image/analysis.py` | Deterministic single-variable layer implemented and tested | Operating-point objective, detuning sweep, probe-power sweep, exposure-time sweep, and sweep summary helper are available for information-versus-destruction analysis. |

## Notebook Section Status

| Section | Status after Version 1 core migration |
|---|---|
| 00 Imports | Notebook-local runtime and plotting setup remain. |
| 01 Parameters | Notebook-local parameter definitions remain authoritative. |
| 02 Atomic Model | Core Thomas-Fermi helper support exists; notebook section remains authoritative. |
| 03 Light-Atom Interaction | Core helper support exists; notebook section remains authoritative. |
| 04 PCI | Core imaging helpers exist; SNR analysis and operating maps remain notebook-local. |
| 05 DGI | Core imaging and camera helpers exist; demonstrations remain notebook-local. |
| 06 Faraday | Core Faraday helper exists; dual-port frame sequence and narrative figures remain notebook-local. |
| 07 Camera | Deterministic and stochastic camera helpers exist; demonstrations remain notebook-local. |
| 08 Shot Noise | Analytical SNR and shot-budget analysis remain notebook-local. |
| 09 Multi-shot Simulation | Deterministic sequence core exists; noisy frame rendering and filmstrips remain notebook-local. |
| 10 Analysis | Deterministic single-variable Faraday optimisation helpers exist; broader narrative analysis, figures, plotting, multi-parameter sweeps, and report-style walkthrough remain notebook-local. |

## Regression Baselines

Current baselines:

```text
regression/baseline/notebook_outputs.json
regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

The `.npz` baselines are tracked in the repository for regression tests but are
excluded from the generated portable code bundle by default.

## Current Validation State

Current validation commands:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

Current results:

```text
pytest -q: 66 passed
notebook section validation: passed
```

## Remaining Notebook-Local Work

The following remain notebook-local and are not part of the closed Version 1
migrated core:

- noisy frame rendering and filmstrips;
- Faraday dual-port frame sequence;
- two-dimensional / three-dimensional optimisation sweeps;
- stochastic noise averaging for optimisation;
- operating maps;
- plotting and figure generation;
- automated optimisation logic;
- experimental RAI calibration;
- broader analysis narrative in `notebook_sections/10_analysis.py`.

## Version 1 Closure Assessment

The Version 1 migrated simulator core is closed for MSc project/report use.

This means the central reusable layers are helperized, exported, and
regression-tested while the original notebook remains the scientific authority.
It does not mean every notebook cell has been migrated or that optimisation and
presentation workflows have been refactored out of the notebook.
