# Optimisation Readiness: Continuous Faraday Imaging

## Purpose

This document outlines what is needed before optimising continuous Faraday
imaging with the Version 1 migrated simulator core.

It is a design note only. It does not implement optimisation, change simulator
code, change helper APIs, modify notebook sections, or alter baselines.

## Optimisation Goal

The optimisation goal is to maximise useful imaging information while
minimising destructiveness.

For this project, "useful information" should stay tied to measurable imaging
quality, such as Faraday signal, per-frame SNR, accumulated SNR, and the number
of frames that remain useful before the condensate is significantly degraded.

"Destructiveness" should remain tied to physically meaningful cost measures:
scattered photons per atom, condensate loss, heating, reabsorption, and the
finite photon budget.

## Candidate Optimisation Variables

Initial optimisation variables could include:

- probe detuning;
- probe intensity or probe power;
- exposure time;
- number of frames;
- imaging mode;
- analyser setting or Faraday readout choice.

For a conservative first optimisation milestone, the variables should be kept
small. A single operating point and a single Faraday readout convention should
be used before attempting a full parameter scan.

## Physical Constraints

Any optimisation should respect:

- scattered photons per atom;
- condensate loss fraction;
- heating of the trapped cloud;
- reabsorption or residual absorption;
- camera shot noise and read noise;
- finite photon budget;
- finite useful frame count;
- the current phenomenological Faraday convention.

The current notebook convention remains:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

`kappa_F` should not be recalibrated inside a first optimisation milestone.
Faraday calibration belongs in a later calibration-specific task.

## Possible Metrics

Useful candidate metrics include:

- SNR per frame;
- RMS accumulated SNR;
- information per scattered photon;
- maximum useful frame number;
- condensate survival fraction;
- Faraday signal-to-destruction ratio;
- accumulated information before a specified loss threshold;
- minimum destructiveness required to reach a target SNR.

The most direct Version 1 metric is a Faraday signal-to-destruction ratio using
deterministic per-frame signal and the existing scattered-photon / multi-shot
bookkeeping.

## Existing Helpers That Support Optimisation

The Version 1 migrated core already provides most lower-level pieces needed for
a controlled optimisation layer:

- `faraday_rotation_angle(...)`
  - Computes the current phenomenological Faraday rotation scalar.

- `scattered_photons_per_atom(...)`
  - Computes scattering cost per atom per probe pulse.

- `reabsorption_fraction(...)`
  - Supports heating / destructiveness estimates from residual absorption.

- `simulate_faraday_image(...)`
  - Produces notebook-equivalent Faraday dark-field and dual-port image arrays
    from a supplied `theta_F` map and pupil.

- `simulate_camera_image(...)`
  - Provides deterministic camera binning and normalisation.

- `simulate_noisy_camera_image(...)`
  - Provides stochastic camera noise with explicit RNG handling when a
    noise-averaged objective is needed.

- `simulate_multishot_sequence(...)`
  - Provides deterministic heating / clean-loss sequence bookkeeping.

- `accumulate_snr(...)`
  - Preserves the notebook RMS accumulated-SNR convention.

These helpers should be composed by an optimisation layer rather than changed.

## Missing Pieces Before Optimisation

The following pieces should be added before real optimisation work:

- an explicit parameter-sweep wrapper;
- a deterministic objective function for Faraday imaging;
- a noise-averaged objective option using explicit RNG seeds;
- a result table or structured result object;
- a plotting layer for operating maps;
- calibration inputs if experimental RAI data becomes available;
- clear selection of which Faraday observable is optimised: dark-field,
  dual-port signal, or another explicitly defined readout.

The first objective should be deterministic. Stochastic or noise-averaged
optimisation can come later once the deterministic objective is validated.

## Conservative Next Milestone

The recommended next milestone is:

```text
Create a deterministic Faraday optimisation objective for one operating point.
```

That milestone should not be a full parameter scan. It should take one chosen
set of operating parameters and return a small, auditable set of values such as:

- Faraday signal scale;
- scattered photons per atom;
- per-frame deterministic SNR estimate;
- accumulated SNR;
- condensate survival fraction;
- signal-to-destruction ratio.

Only after this objective is checked against notebook expectations should the
project add parameter sweeps, plotting, stochastic averaging, or optimisation
search algorithms.

## Boundary Conditions

Near-term optimisation work should not:

- redesign the Faraday model;
- change `kappa_F`;
- alter helper APIs;
- modify notebook physics;
- mix plotting into the objective function;
- combine RAI calibration and optimisation in the same milestone;
- implement beyond-Thomas-Fermi state models inside the optimisation layer.

The optimisation layer should sit above the Version 1 migrated core and consume
its helpers as stable building blocks.
