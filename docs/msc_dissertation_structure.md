# MSc Dissertation Structure Draft

## 1. Proposed Title

**Simulation and Optimisation of Continuous Non-Destructive Imaging for an
Ultracold `166Er` Bose-Einstein Condensate**

Alternative shorter title:

**Continuous Non-Destructive Imaging of an Ultracold `166Er` Condensate:
Simulation, Validation, and Faraday Optimisation**

## 2. Abstract Outline

The abstract should cover:

- the scientific motivation for continuous non-destructive imaging of
  ultracold gases;
- the need to balance useful imaging information against destructiveness;
- the use of the original notebook as the authoritative scientific
  implementation;
- the conservative migration of notebook-equivalent behaviour into a tested
  helper package;
- implementation of Atomic Model, Light-Atom Interaction, Imaging, Camera,
  Stochastic Camera Noise, Deterministic Multi-shot, and deterministic Faraday
  optimisation layers;
- validation through regression baselines, tests, and notebook-section
  synchronisation checks;
- deterministic Faraday operating-point analysis using detuning, probe
  intensity, and exposure-time sweeps;
- limitations and future directions, including calibration with experimental
  RAI data and beyond-Thomas-Fermi states.

## 3. Chapter Structure

### Chapter 1: Introduction

- Motivation for non-destructive imaging in ultracold-atom experiments.
- Why continuous imaging is useful for dynamic condensates.
- Main trade-off: information gained versus destructiveness.
- Project scope and constraints.
- Summary of contributions.

### Chapter 2: Physical Background

- Bose-Einstein condensates and the Thomas-Fermi approximation.
- Relevant properties of ultracold `166Er`.
- Light-atom interaction concepts used in the simulator.
- PCI, DGI, and Faraday imaging at a conceptual level.
- Scattering, heating, reabsorption, and condensate loss.

### Chapter 3: Original Notebook Simulator

- Role of the original notebook as the authoritative implementation.
- Overview of notebook sections and scientific workflow.
- Existing notebook conventions and assumptions.
- Current Faraday convention:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

- Limitations of a notebook-only implementation for validation and extension.

### Chapter 4: Version 1 Migrated Simulator Core

- Conservative migration strategy.
- Helper package architecture.
- Layer-by-layer description:
  Atomic Model, Light-Atom Interaction, Imaging, Camera, Stochastic Camera
  Noise, Deterministic Multi-shot Core.
- Public helper API and separation of concerns.
- Confirmation that the helper package preserves notebook-equivalent behaviour
  and does not replace the original notebook as scientific authority.

### Chapter 5: Validation and Reproducibility

- Regression baseline strategy.
- Notebook-section validation.
- Test suite structure.
- Code bundle and reproducibility policy.
- Large-file handling and baseline management.
- Current validation status.

### Chapter 6: Deterministic Faraday Optimisation Layer

- Motivation for optimisation above the migrated core.
- Operating-point objective.
- One-dimensional detuning, probe-intensity, and exposure-time sweeps.
- Sweep summary helper.
- Interpretation of metrics such as `signal_per_scattered_photon` and
  `information_per_scattered_photon`.
- Limitations of the current deterministic layer.

### Chapter 7: Results

- Demonstration of migrated helper behaviour.
- Representative PCI/DGI and Faraday imaging baselines.
- Camera and multi-shot deterministic outputs.
- Faraday optimisation examples.
- Information-versus-destruction comparisons.
- Discussion of what the current results do and do not prove.

### Chapter 8: Discussion

- Scientific interpretation of the simulator outputs.
- Strengths of the migrated architecture.
- Limits of the Thomas-Fermi and phenomenological Faraday assumptions.
- Impact of fixed `kappa_F`.
- Practical implications for experimental imaging choices.
- Risks in interpreting deterministic optimisation without calibration.

### Chapter 9: Future Work

- Experimental RAI calibration.
- Droplets, supersolids, and mixtures.
- `kappa_F` calibration.
- Stochastic optimisation.
- Two-dimensional and three-dimensional sweeps.
- Report-ready plotting and operating maps.

### Chapter 10: Conclusion

- Summary of simulator migration.
- Summary of validation and reproducibility contributions.
- Summary of Faraday optimisation readiness.
- How the project supports future experimental and theoretical work.

## 4. Key Research Question

The central research question can be framed as:

**How can continuous non-destructive imaging of an ultracold `166Er`
Bose-Einstein condensate be simulated and optimised to maximise useful
information while minimising destructiveness?**

Supporting questions:

- How much useful imaging signal is obtained per scattered photon?
- How do detuning, probe intensity, and exposure time affect the
  information-versus-destruction trade-off?
- How can a notebook-based simulator be made reproducible and maintainable
  without changing its scientific behaviour?
- Which assumptions must be calibrated or extended before experimental
  optimisation claims can be made?

## 5. Methodology Section Outline

The methodology chapter should explain:

- the original notebook-first development process;
- the decision to preserve notebook-equivalent behaviour;
- the staged migration strategy;
- how helper functions were extracted by layer;
- how frozen APIs prevented accidental physics changes;
- how regression baselines were generated and used;
- how tests were added for each migrated layer;
- how Faraday optimisation helpers were added only after the imaging, camera,
  and multi-shot core was stable.

Suggested subsections:

- Notebook audit and migration boundary.
- Helper extraction strategy.
- Baseline generation.
- Regression testing.
- Deterministic optimisation workflow.
- Documentation and reproducibility workflow.

## 6. Simulator Architecture Section Outline

The architecture section should describe the current flow:

```text
Original notebook, authoritative scientific implementation
  -> notebook section exports
  -> migrated helper package
  -> regression tests and baselines
  -> deterministic optimisation layer
```

Layer stack:

```text
Atomic Model
  -> Light-Atom Interaction
  -> Imaging
  -> Camera / Stochastic Camera Noise
  -> Deterministic Multi-shot Core
  -> Deterministic Faraday Optimisation
  -> Notebook-local plotting, narrative analysis, and future extensions
```

Key points:

- lower layers remain explicit and reusable;
- helper inputs avoid hidden notebook globals;
- the original notebook remains the scientific source of truth;
- the helper package preserves notebook-equivalent behaviour;
- optimisation sits above the migrated simulator core rather than changing it.

## 7. Physics Model Section Outline

The physics model section should include:

- Thomas-Fermi condensate model;
- recoil and atomic quantities;
- scalar dispersive phase shift;
- residual optical depth;
- scattered photons per atom;
- reabsorption fraction;
- PCI and DGI scalar imaging conventions;
- Fourier propagation and pupil conventions;
- Faraday phenomenological rotation model;
- camera binning, normalisation, and stochastic noise;
- deterministic multi-shot heating / clean-loss bookkeeping;
- accumulated SNR convention.

Important boundary statement:

The current Faraday model remains phenomenological. It uses:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

The project does not yet introduce a microscopic circular-transition Faraday
model or recalibrate `kappa_F`.

## 8. Validation / Reproducibility Section Outline

This section should describe:

- unit and regression tests;
- notebook-section validation script;
- stored notebook-output baseline;
- PCI/DGI imaging baseline;
- Faraday imaging baseline;
- deterministic tests for imaging, camera, multi-shot, and optimisation
  helpers;
- reproducible code bundle generation;
- large-file hygiene policy.

Current validation state:

```text
pytest -q: 66 passed
notebook section validation: passed
```

The section should emphasise that validation checks notebook-equivalent helper
behaviour. It does not claim that the underlying physical model is complete or
experimentally calibrated.

## 9. Faraday Optimisation Section Outline

This section should explain:

- why Faraday optimisation is framed as information versus destructiveness;
- why maximum image brightness is not the correct objective by itself;
- how `evaluate_faraday_operating_point(...)` estimates deterministic signal,
  scattering, destructiveness, and information-efficiency quantities;
- how detuning sweeps compare dispersive signal against scattering cost;
- how probe-intensity sweeps compare signal efficiency against destructive
  scaling;
- how exposure-time sweeps compare acquisition time against scattering cost;
- how `summarise_faraday_sweep(...)` selects the best parameter according to an
  explicit metric.

Metrics to discuss:

- `signal_per_scattered_photon`;
- `information_per_scattered_photon`;
- `signal_to_destruction`;
- `estimated_per_frame_snr`.

Limitations:

- deterministic only;
- single-variable sweeps only;
- no plotting layer yet;
- no stochastic averaging yet;
- no RAI calibration yet;
- fixed phenomenological `kappa_F`.

## 10. Results Section Plan

Possible result figures and tables:

- architecture diagram of migrated simulator layers;
- table of migrated helpers and validation status;
- representative Thomas-Fermi density or column-density profile;
- PCI/DGI baseline image arrays or summary statistics;
- Faraday baseline dark-field and dual-port image outputs;
- camera pipeline demonstration;
- deterministic multi-shot sequence plot or table;
- Faraday detuning sweep table;
- Faraday probe-intensity sweep table;
- Faraday exposure-time sweep table;
- summary table comparing best operating points under selected metrics.

Result narrative:

- first demonstrate that the simulator core is validated and reproducible;
- then show representative imaging behaviour;
- then present deterministic Faraday optimisation as an initial
  information-versus-destruction analysis;
- clearly separate validated simulation outputs from future experimental
  calibration claims.

## 11. Discussion and Limitations

Discussion topics:

- benefits of moving from notebook-only code to a tested helper architecture;
- how layer separation supports reproducibility and future optimisation;
- interpretation of information per scattered photon;
- practical meaning of deterministic sweep optima;
- why the best information-efficiency point may not be the brightest image;
- impact of camera noise and multi-shot destructiveness on future optimisation.

Limitations:

- original notebook remains authoritative, so helper package is not an
  independent new model;
- Thomas-Fermi state model does not yet cover droplets, supersolids, or
  mixtures;
- Faraday model is phenomenological;
- `kappa_F` remains fixed at `1.0`;
- no experimental RAI calibration yet;
- no stochastic optimisation or noise averaging yet;
- no two-dimensional or three-dimensional parameter sweeps yet;
- plotting and final operating maps remain future presentation work.

## 12. Future Work

### RAI Calibration

Use experimental RAI data to calibrate atom number, cloud widths, optical depth
scale, magnification, effective pixel size, camera response, detuning offset,
and residual absorption. Calibration parameters should be stored in explicit
external config files, not hard-coded into lower-level helpers.

### Droplets, Supersolids, and Mixtures

Add new state/profile generators rather than modifying
`build_thomas_fermi_state(...)`. Practical first steps include phenomenological
profiles or external density-map input from experiment or separate simulations.

### `kappa_F` Calibration

Treat `kappa_F` calibration as a separate Faraday-specific calibration
milestone. It should not be changed inside migration or deterministic
optimisation helper work.

### Stochastic Optimisation

Extend deterministic objectives with explicit-RNG noise averaging after the
deterministic optimisation layer is stable. This should use the existing
stochastic camera noise policy and avoid hidden global random state.

### 2D / 3D Sweeps

Design multi-parameter sweeps before implementation. Candidate combinations
include:

- detuning versus probe intensity;
- detuning versus exposure time;
- probe intensity versus exposure time;
- detuning, intensity, and exposure time together.

These should produce structured result tables before plotting or automated
search algorithms are added.

## Closing Note

The dissertation should repeatedly distinguish three levels:

1. The original notebook is the authoritative scientific implementation.
2. The migrated helper package preserves notebook-equivalent behaviour and
   makes the simulator testable and maintainable.
3. The deterministic Faraday optimisation layer provides the first
   optimisation-ready interface for studying information gained versus
   destructiveness.
