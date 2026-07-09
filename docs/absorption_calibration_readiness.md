# Absorption Image Calibration Readiness

## Purpose

This document describes the first conservative calibration-readiness layer for
future laboratory RAI / absorption-image data.

The current implementation adds deterministic preprocessing and observable
extraction helpers only. It does not load lab-specific file formats, fit
parameters, change simulator physics, replace notebook behaviour, calibrate
`kappa_F`, or implement mixture cross-talk correction.

## Current Helpers

The calibration readiness helpers are:

- `compute_optical_density(...)`
- `integrate_optical_density(...)`
- `estimate_cloud_moments(...)`
- `extract_absorption_observables(...)`

`compute_optical_density(...)` computes:

```text
OD = -log((atom_image - dark_image) / (probe_image - dark_image))
```

with safe clipping controlled by `epsilon` to avoid division by zero or
`log(0)`.

`integrate_optical_density(...)` returns the summed OD, optionally scaled by a
pixel area.

`estimate_cloud_moments(...)` estimates:

- peak OD;
- integrated OD;
- centre position;
- RMS width in x;
- RMS width in y.

`extract_absorption_observables(...)` combines OD calculation and moment
extraction into a small dictionary for later calibration workflows.

## Role of Absorption / RAI Data

Experimental RAI or absorption images can provide destructive calibration data
for the non-destructive imaging simulator. They can be used as a ground-truth
reference for quantities such as:

- atom number or atom-number proxy;
- cloud centre;
- cloud widths;
- peak optical density;
- integrated optical density;
- column-density scale;
- magnification;
- effective pixel size;
- camera offset and gain;
- probe and reference intensity balance;
- residual absorption scale.

This calibration should live above the migrated simulator core. It should not
silently replace notebook constants or modify lower-level physics helpers.

## External Calibration Policy

Future fitted calibration parameters should be stored in explicit external
files, for example:

```text
calibration/rai_calibration_example.json
calibration/lab_run_YYYYMMDD.toml
```

The simulator should load those parameters deliberately in a calibration or
analysis layer. Values should not be hard-coded into Atomic Model,
Light-Atom Interaction, Imaging, Camera, Multi-shot, or optimisation helpers.

## Notebook Authority

The original notebook remains the authoritative scientific implementation.

The calibration helpers are a support layer for future data ingestion and
parameter estimation. They do not replace notebook physics, change helper API
behaviour, alter regression baselines, or modify notebook sections.

## Faraday Calibration Boundary

The current Faraday convention remains:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

This calibration-readiness milestone does not calibrate `kappa_F`. Any future
`kappa_F` calibration should be a separate Faraday-specific calibration
milestone using explicit experimental observables and validation data.

## Future Workflow

A conservative future absorption-calibration workflow would be:

1. Load atom, probe, and dark images from a lab-specific file format.
2. Compute OD maps using `compute_optical_density(...)`.
3. Extract observables using `extract_absorption_observables(...)`.
4. Compare observables to simulator predictions.
5. Fit a small set of calibration parameters.
6. Store fitted parameters in external config files.
7. Validate on held-out absorption images.
8. Use calibrated parameters in non-destructive PCI, DGI, and Faraday
   simulations.

## Cross-Talk and Mixture Future Work

Mixture and cross-talk calibration are future extensions. They may require:

- component-resolved OD maps;
- a channel response matrix;
- species-specific or transition-specific cross sections;
- spin- or polarisation-resolved calibration parameters;
- polarisation-channel imbalance correction;
- correction for probe leakage or imperfect analyser extinction;
- validation against held-out absorption images.

These extensions should be additive. They should not change the current
single-component Thomas-Fermi simulator path or the existing notebook-equivalent
helper APIs.

## Current Limitations

The current calibration-readiness layer:

- does not load real lab data;
- does not implement fitting;
- does not infer atom number directly;
- does not calibrate magnification or camera gain by itself;
- does not perform uncertainty estimation;
- does not perform Faraday recalibration;
- does not handle mixture cross-talk.

It provides only deterministic building blocks for future experimental
calibration.
