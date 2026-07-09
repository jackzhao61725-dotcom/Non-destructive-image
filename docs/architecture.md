# Architecture Overview

This document summarises the repository architecture after completion of the
Version 1 migrated simulator core.

It is documentation only. It does not define new physics, rename APIs, or move
code.

## Scientific Authority

The original notebook remains the authoritative scientific implementation:

```text
1 calculations revised 2  multishot  6  extended.ipynb
```

The helper package in `src/non_destructive_image/` is a conservative support
layer. It preserves notebook-equivalent behaviour for migrated core operations
and is regression-tested against stored notebook outputs and numerical imaging
baselines.

Version 1 migration did not introduce a physics redesign.

## Logical Layers

```text
Atomic Model
  |
  v
Light-Atom Interaction
  |
  v
Imaging
  |
  v
Camera / Stochastic Camera Noise
  |
  v
Deterministic Multi-shot Core
  |
  v
Notebook-local Analysis / Plotting / Optimisation
```

The reusable Version 1 core is helperized through the Multi-shot layer. The
analysis, plotting, noisy frame-sequence presentation, and optimisation
workflows remain notebook-local.

## Helper Package Responsibilities

| Module | Public helpers | Responsibility |
|---|---|---|
| `atomic_model.py` | `ThomasFermiState`, `build_thomas_fermi_state`, `recoil_quantities` | Thomas-Fermi state construction and recoil quantities. |
| `profiles.py` | `thomas_fermi_profile_2d` | Reusable 2D Thomas-Fermi profile expression. |
| `light_atom.py` | `dimensionless_detuning`, `scalar_phase_shift`, `residual_optical_depth`, `intensity_at_atoms`, `scattered_photons_per_atom`, `faraday_rotation_angle`, `reabsorption_fraction` | Optical interaction quantities derived from notebook-equivalent formulas. |
| `fourier.py` | `propagate_scattered_field` | Notebook FFT/pupil propagation operation. |
| `imaging.py` | `simulate_fourier_image`, `simulate_pci_image`, `simulate_dgi_image`, `simulate_faraday_image` | Shared Fourier imaging core and PCI/DGI/Faraday orchestration helpers. |
| `camera.py` | `bin_to_camera_pixels`, `add_camera_noise`, `normalize_camera_counts`, `simulate_camera_image`, `simulate_noisy_camera_image` | Deterministic camera binning/normalisation and stochastic camera-noise orchestration with explicit RNG handling. |
| `multishot.py` | `simulate_multishot_sequence`, `accumulate_snr` | Deterministic multi-shot sequence bookkeeping, heating / clean-loss updates, and RMS accumulated-SNR convention. |

## Current Data Flow

```text
Physical parameters and notebook constants
  |
  v
Atomic state and Thomas-Fermi profile
  |
  v
Phase, scattering, Faraday rotation, and reabsorption quantities
  |
  v
PCI / DGI / Faraday image formation
  |
  v
Deterministic or stochastic camera image
  |
  v
Deterministic multi-shot sequence bookkeeping
```

The migrated helper stack avoids notebook globals. Inputs are explicit, and
future calibration or optimisation layers should call the helpers rather than
changing lower-level formulas.

## Version 1 Closure

The Version 1 migrated core is closed for MSc project/report use:

- Atomic Model helper layer is implemented and tested.
- Light-Atom Interaction helper layer is implemented and tested.
- Shared Fourier imaging core is implemented and tested.
- PCI orchestration helper is implemented and tested.
- DGI orchestration helper is implemented and tested.
- Faraday orchestration helper is implemented and tested.
- Deterministic camera pipeline is implemented and tested.
- Stochastic camera noise helper is implemented and tested with explicit RNG.
- Deterministic multi-shot sequence core is implemented and tested.

Current validation status:

```text
pytest -q: 37 passed
notebook section validation: passed
```

## Remaining Notebook-Local Work

The following workflows are intentionally outside the closed Version 1 core:

- noisy frame rendering and filmstrips;
- Faraday dual-port frame sequence;
- detuning sweep and operating maps;
- plotting and figure generation;
- optimisation logic;
- broader narrative analysis in `notebook_sections/10_analysis.py`.

These can be migrated later, but they should not be mixed into the current
Version 1 core closure.

## Architecture Rules For Future Work

1. Keep the original notebook as the scientific authority unless a future
   milestone explicitly replaces a section after validation.
2. Do not rename or change frozen public helpers during report integration.
3. Keep `kappa_F` explicit in Faraday-related work until a separate
   calibration milestone changes the model.
4. Keep calibration and optimisation layers above the migrated core.
5. Add new state models, RAI calibration, and beyond-Thomas-Fermi extensions
   additively rather than replacing the Thomas-Fermi Version 1 path.

## Related Documents

- `docs/migration_status.md`
- `docs/version_1_migrated_core_summary.md`
- `docs/extension_roadmap.md`
