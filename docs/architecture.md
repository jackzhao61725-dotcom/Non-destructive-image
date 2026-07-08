# Architecture Overview

This document reviews the current repository architecture after the Atomic Model and Light-Atom Interaction migrations. It is an architecture-freeze review only: it does not define new physics, rename APIs, or move code.

## 1. Current repository architecture

```text
Original Notebook (authoritative scientific reference)
  |
  | exported and audited
  v
Notebook Section Exports (migration surface)
  |
  | selected formulas replaced by explicit helper calls
  v
Helper Package (tested formula layer)
  |
  | compared by regression tests / future numerical baselines
  v
Future Refactored Notebook or Scripts (not implemented yet)
```

The original notebook, `1 calculations revised 2  multishot  6  extended.ipynb`, remains the scientific authority. The section exports are migration targets, and the helper package is an explicit, testable implementation surface for formulas that have already been extracted.

## 2. Current logical layers

```text
00 Imports / runtime setup
  |
01 Parameters and constants
  |
02 Atomic Model
  |  - trap frequencies
  |  - Thomas-Fermi state
  |  - density, radii, column densities
  v
03 Light-Atom Interaction
  |  - detuning
  |  - scalar phase
  |  - residual optical depth
  |  - scattering per pulse
  |  - phenomenological Faraday rotation angle
  v
04-06 Imaging Models
  |  - PCI
  |  - DGI
  |  - Faraday field/Jones propagation
  v
07 Camera
  |  - camera binning
  |  - digitisation / normalisation conventions
  v
08 Shot Noise and SNR
  |
09 Multi-shot Simulation
  |  - heating/loss sequence
  |  - frame-by-frame state update
  v
10 Analysis / plotting / comparisons
```

This layering is appropriate for the remaining migration work. The stable boundary should be:

```text
Atomic Model -> Light-Atom Interaction -> Imaging -> Camera/Noise -> Multi-shot/Analysis
```

A separate Fourier/propagation layer may sit underneath PCI and DGI because PCI/DGI share a Fourier-plane propagation pattern, while the current Faraday model uses Jones/polarisation propagation rather than Fourier-plane masks.

## 3. Helper package responsibilities

| Module | Current public helpers | Responsibility | Boundary assessment |
|---|---|---|---|
| `atomic_model.py` | `ThomasFermiState`, `build_thomas_fermi_state`, `recoil_quantities` | Derived Thomas-Fermi state and recoil quantities from explicit physical inputs. | Correct module. Atomic quantities are separated from optical calculations. |
| `profiles.py` | `thomas_fermi_profile_2d` | 2D Thomas-Fermi column-profile shape used by imaging demonstrations. | Correct for now. Could later be folded into an imaging/profile submodule, but no change required before freeze. |
| `light_atom.py` | `dimensionless_detuning`, `scalar_phase_shift`, `residual_optical_depth`, `intensity_at_atoms`, `scattered_photons_per_atom`, `faraday_rotation_angle`, `reabsorption_fraction` | Converts column density, detuning, optical constants, and probe settings into optical interaction quantities. | Correct module. `reabsorption_fraction` is borderline because it feeds heating/loss, but it is still derived from residual OD and belongs here for now. |
| `fourier.py` | `propagate_scattered_field` | Shared FFT/pupil operation for Fourier-plane imaging paths. | Correct module. It supports PCI/DGI Fourier-core helpers without coupling to atom physics. |
| `imaging.py` | `simulate_fourier_image` | Shared coherent PCI/DGI Fourier-imaging core: object field, pupil/mask propagation, reference-field recombination, optional intensity. | Correct module. It is intentionally below PCI/DGI-specific orchestration and above the raw FFT helper. |
| `camera.py` | `bin_to_camera_pixels`, `add_camera_noise`, `normalize_camera_counts` | Camera binning/noise/count normalisation patterns. | Correct module. It should remain independent of physical image generation. |

## 4. Public helper API review

### Atomic Model

- `build_thomas_fermi_state(...)`
  - Naming: clear and descriptive.
  - Scope: returns the derived Thomas-Fermi state while keeping inputs explicit.
  - Responsibility: correct; it does not include optical or imaging calculations.
  - Freeze assessment: stable enough for future migration.

- `ThomasFermiState`
  - Naming: clear.
  - Scope: appropriate container for derived atomic quantities.
  - Responsibility: correct; it gives downstream code a structured state while preserving scalar fields.
  - Freeze assessment: stable enough, although future baselines should confirm every field numerically.

- `recoil_quantities(...)`
  - Naming: clear.
  - Scope: returns recoil energy, recoil temperature, and recoil velocity.
  - Responsibility: acceptable in `atomic_model.py` because recoil is atom/transition dependent and feeds later heating checks.
  - Freeze assessment: stable enough.

### Thomas-Fermi profile

- `thomas_fermi_profile_2d(...)`
  - Naming: clear.
  - Scope: returns only the profile expression, not a full image.
  - Responsibility: acceptable as a profile helper used by imaging code.
  - Freeze assessment: stable enough, but future imaging migration should decide whether profile generation belongs under imaging or remains as a shared profile utility.

### Light-Atom Interaction

- `dimensionless_detuning(...)`
  - Naming: clear.
  - Scope: one formula, no hidden globals.
  - Responsibility: correct.
  - Freeze assessment: stable.

- `scalar_phase_shift(...)`
  - Naming: clear and physically meaningful.
  - Scope: converts column density and detuning into scalar phase.
  - Responsibility: correct; no imaging assumptions.
  - Freeze assessment: stable.

- `residual_optical_depth(...)`
  - Naming: clear.
  - Scope: residual OD at detuning.
  - Responsibility: correct.
  - Freeze assessment: stable.

- `intensity_at_atoms(...)`
  - Naming: clear.
  - Scope: converts probe power and beam diameter into intensity using the notebook convention.
  - Responsibility: acceptable in Light-Atom Interaction because it is needed before scattering.
  - Freeze assessment: stable.

- `scattered_photons_per_atom(...)`
  - Naming: clear.
  - Scope: scattering per atom per pulse.
  - Responsibility: correct.
  - Freeze assessment: stable.

- `faraday_rotation_angle(...)`
  - Naming: clear.
  - Scope: current phenomenological scalar rotation angle only.
  - Responsibility: correct for the documented notebook model.
  - Freeze assessment: stable only if reviewers accept that `kappa_f` remains an explicit calibration placeholder.

- `reabsorption_fraction(...)`
  - Naming: clear.
  - Scope: computes reabsorption from residual OD along principal axes.
  - Responsibility: acceptable but borderline; it supports heating/loss logic while depending on Light-Atom OD.
  - Freeze assessment: acceptable for now; future destruction/heating modules may wrap it rather than move it.

### Fourier propagation

- `propagate_scattered_field(...)`
  - Naming: clear.
  - Scope: one FFT/pupil/IFFT operation.
  - Responsibility: correct; it does not know about PCI, DGI, or atoms.
  - Freeze assessment: stable enough for future PCI/DGI migration.

### Camera

- `bin_to_camera_pixels(...)`
  - Naming: clear.
  - Scope: camera binning only.
  - Responsibility: correct.
  - Freeze assessment: stable.

- `add_camera_noise(...)`
  - Naming: clear.
  - Scope: stochastic camera noise recipe.
  - Responsibility: correct, but future baselines must specify seeds carefully.
  - Freeze assessment: stable if random-seed policy remains explicit.

- `normalize_camera_counts(...)`
  - Naming: clear.
  - Scope: count normalisation only.
  - Responsibility: correct.
  - Freeze assessment: stable.

## 5. Data flow review

Current intended scientific data flow:

```text
Parameters / constants
  |
  v
Atomic Model
  - atom number
  - scattering length
  - trap frequencies
  - mass
  |
  v
Thomas-Fermi State
  - chemical potential
  - peak density
  - radii
  - column density
  |
  v
Light-Atom Interaction
  - detuning
  - scalar phase
  - residual optical depth
  - scattering per atom per shot
  - Faraday rotation angle
  |
  v
Imaging Model
  - PCI object field and Fourier-plane phase dot
  - DGI object field and dark mask
  - Faraday polarisation/Jones propagation
  |
  v
Camera / Noise
  - camera binning
  - photon/read noise
  - normalised image
  |
  v
Multi-shot Simulation
  - heating/loss update
  - repeated frame generation
  |
  v
Analysis
  - SNR
  - plots
  - comparison tables
```

The flow is scientifically clean at the conceptual level. The remaining implementation coupling is caused by notebook globals and by plotting/analysis cells mixing calculations with display code.

## 6. Coupling and boundary risks

Current unnecessary coupling:

- Notebook sections still depend on global names and execution order.
- `03_light_atom_interaction.py` keeps wrapper functions that bind helper inputs to notebook globals.
- `06_faraday.py` still mixes scalar Faraday rotation tables with Faraday field propagation and plotting.
- `09_multishot_simulation.py` still has direct phase/OD calculations and state updates inside one notebook-style section.
- Plotting and analysis code are still interleaved with scientific calculations in several exported sections.

These couplings are acceptable before Imaging migration as long as future migrations remain gated and module-scoped.

## 7. Future extensibility assessment

### PCI

The current API can support PCI migration. PCI needs scalar phase maps and Fourier propagation. `scalar_phase_shift(...)`, `thomas_fermi_profile_2d(...)`, and `propagate_scattered_field(...)` provide the necessary lower-level ingredients, but PCI-specific assembly should become its own future imaging helper.

### DGI

The current API can support DGI migration. DGI shares scalar phase and Fourier propagation with PCI, but has a different Fourier-plane mask. No public API redesign is required before DGI migration.

### Faraday

The current API can support the notebook's current Faraday model because `faraday_rotation_angle(...)` exposes `kappa_f` explicitly and does not hide the calibration placeholder. It does not solve the theoretical calibration of `kappa_f`; that remains a physics-review issue, not an architecture blocker.

### Optimisation over detuning, intensity, and exposure time

The current API is acceptable for future optimisation because detuning, power/intensity, pulse duration, optical constants, and calibration factors are explicit inputs. However, optimisation may become verbose if every helper receives many scalar arguments. A future ergonomics-only layer may introduce parameter objects, but the current public formula helpers should remain stable.

## 8. Recommended stable architecture

```text
src/non_destructive_image/

  atomic_model.py
    - atomic state and recoil quantities

  profiles.py
    - reusable Thomas-Fermi profile shapes

  light_atom.py
    - detuning
    - scalar phase
    - optical depth
    - scattering
    - phenomenological Faraday rotation
    - reabsorption from residual OD

  fourier.py
    - generic Fourier/pupil propagation operation

  camera.py
    - camera binning
    - camera noise
    - count normalisation

  imaging.py
    - shared PCI/DGI coherent Fourier-imaging core

  future imaging module(s) [not yet implemented]
    - PCI-specific image construction
    - DGI-specific image construction
    - Faraday image construction

  future simulation module(s) [not yet implemented]
    - multi-shot update loop
    - heating/loss sequence
    - optimisation workflows
```

## 9. Architecture freeze recommendation

**APPROVE ARCHITECTURE**, with constraints.

The current helper API is sufficiently clean and stable to support the next migration stage, provided future work follows these rules:

1. Do not rename current public helpers during Imaging migration.
2. Keep physics formulas in their current modules unless a later milestone explicitly approves a module-boundary change.
3. Keep `kappa_f` explicit in Faraday helpers until a separate physics-calibration milestone changes the model.
4. Add future imaging helpers above the current formula/propagation helpers rather than changing the lower-level helpers.
5. Continue preserving notebook wrapper names during migration until numerical baselines prove equivalence.

## 10. Self review

Is the current architecture mature enough that future migration can proceed without changing public APIs?

**YES**, with limitations.

The public API is stable enough for PCI/DGI/Faraday imaging migration. The known limitations are not blockers:

- A future parameter-object layer may improve ergonomics, but it is not required for correctness.
- `reabsorption_fraction(...)` sits near a module boundary, but moving it now would be churn rather than a scientific improvement.
- Full numerical baselines are still missing in this environment, so architecture freeze should not be interpreted as full scientific validation.

## 11. Future Extensions

Future experimental calibration and beyond-Thomas-Fermi state support should be
added after the notebook-equivalent migration is stable. The current Version 1
helpers should remain focused on preserving notebook behaviour.

See `docs/extension_roadmap.md` for the proposed roadmap for RAI-data-based
calibration, droplets, supersolids, mixtures, and later calibrated
optimisation.
