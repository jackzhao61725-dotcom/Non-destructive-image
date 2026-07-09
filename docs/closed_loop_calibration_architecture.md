# Closed-Loop Calibration Architecture

## Purpose

The simulator is not intended to be a one-way prediction tool. It is intended
as a calibration-aware optimisation framework for continuous non-destructive
imaging.

The project is a closed-loop framework: a notebook-equivalent physical model
informs experimental design decisions; the Faraday imaging system generates
measurements; destructive absorption images provide calibration feedback; and
`kappa_F` is ultimately fixed by experiment. The simulator is therefore not a
one-way predictor, but a calibration-aware optimisation tool.

This document is an architecture note only. It does not implement fitting,
change simulator physics, change helper APIs, modify notebook sections, update
baselines, or recalibrate `kappa_F`.

## Closed-Loop Architecture

The intended loop is:

```text
notebook-equivalent physical model
  -> operating-point design decisions
  -> Faraday imaging experiment
  -> destructive absorption / RAI calibration feedback
  -> density, OD, camera, and kappa_F calibration
  -> updated simulation and optimisation
```

The design decision stage may choose:

- detuning;
- probe intensity or probe power;
- exposure time;
- frame number;
- analyser or readout convention in later work.

The experiment stage uses the Faraday imaging system to collect
non-destructive measurements. Destructive absorption / RAI images then provide
ground-truth feedback for calibration. The calibrated simulator can be reused
to update operating-point choices and improve the next experiment.

## Absorption / RAI Observables

Absorption and RAI images can provide stable destructive reference data. The
first useful observables include:

- optical density;
- integrated optical density;
- cloud centre;
- cloud widths;
- peak optical density.

These observables can support calibration of:

- density scale;
- Thomas-Fermi radii or empirical cloud size;
- optical-depth scale;
- magnification and effective pixel size;
- camera offset and gain;
- residual absorption scale.

The calibrated density and OD parameters can then improve PCI, DGI, and
Faraday simulations without changing the notebook-equivalent lower-level
physics.

## Current Absorption Calibration Readiness Layer

The repository already contains first-stage absorption calibration helpers:

- `compute_optical_density(...)`
- `integrate_optical_density(...)`
- `estimate_cloud_moments(...)`
- `extract_absorption_observables(...)`

These helpers support deterministic preprocessing and observable extraction.
They do not yet implement:

- real-data file loaders;
- parameter fitting;
- uncertainty estimation;
- Faraday `kappa_F` fitting;
- mixture cross-talk correction.

Their role is to prepare a clean calibration input layer that can later feed a
fitting pipeline.

## `kappa_F` Interpretation

The Version 1 notebook-equivalent Faraday model uses:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

`kappa_F = 1.0` is a Version 1 placeholder convention. It should not be changed
during notebook-equivalent migration, regression baseline generation, or the
current deterministic optimisation layer.

In the experimental calibration stage, `kappa_F` becomes a fitted or calibrated
parameter. This is not a contradiction. It is the intended transition from an
idealised notebook-equivalent simulator to an experimentally calibrated
simulator.

The correct boundary is:

- keep `kappa_F = 1.0` fixed in Version 1 migration and deterministic helper
  validation;
- calibrate `kappa_F` later using explicit experimental observables and
  held-out validation data;
- store calibrated values in external calibration files rather than silently
  changing lower-level helper defaults.

## Relation To Faraday Optimisation

The current deterministic Faraday optimisation layer can compare operating
points using metrics such as:

- `signal_per_scattered_photon`;
- `information_per_scattered_photon`;
- `signal_to_destruction`;
- `estimated_per_frame_snr`.

Those metrics are useful before calibration, but their experimental meaning
improves after absorption / RAI data fixes density scale, OD scale, camera
response, and eventually `kappa_F`.

The intended workflow is therefore:

1. Use the notebook-equivalent model to choose candidate operating points.
2. Run the experiment and collect Faraday measurements.
3. Collect destructive absorption / RAI reference images.
4. Extract OD and cloud observables.
5. Fit calibration parameters.
6. Update the simulator inputs.
7. Repeat optimisation with calibrated parameters.

## Relation To Cross-Talk And Mixtures

The current Version 1 path assumes a simpler single-component calibration path.
Future mixture systems may require component-resolved calibration and explicit
cross-talk correction.

Cross-talk belongs in the calibration layer, not in the notebook-equivalent
simulator core.

A possible future model is:

```text
measured channel vector = response matrix * component density vector + noise
```

Future correction may require:

- component-resolved OD maps;
- channel response matrices;
- species-specific or transition-specific cross sections;
- spin- or polarisation-resolved calibration parameters;
- polarisation-channel imbalance correction;
- validation against held-out absorption images.

These extensions should be additive. They should not replace the current
single-component Thomas-Fermi notebook-equivalent path.

## MSc Relevance

This architecture makes the project more than a code refactor and more than a
forward image simulator.

It frames the simulator as a calibration-aware optimisation framework for
continuous non-destructive imaging. The model proposes operating points, the
experiment tests them, destructive absorption images provide calibration
feedback, and the calibrated simulator improves subsequent design decisions.

This is directly connected to the core MSc question: how to maximise useful
imaging information while minimising destructiveness under experimental
constraints.

## Current Boundary

The repository currently does not include:

- a real-data absorption / RAI loader;
- a full fitting pipeline;
- `kappa_F` fitting;
- cross-talk correction;
- mixture response matrices;
- calibrated Faraday optimisation;
- changes to notebook-equivalent physics.

The current implementation provides only the validated simulator core,
deterministic single-variable Faraday optimisation helpers, and first-stage
absorption observable extraction helpers needed for future calibration work.
