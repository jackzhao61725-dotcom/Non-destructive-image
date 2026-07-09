# Notebook-Aligned Recovery Status

## Current Recovery Path

The canonical notebook-aligned recovery path now covers:

```text
parameters -> Thomas-Fermi condensate -> projected profile -> column-density map -> scalar phase map
```

This remains a recovery of the historical Version 1 notebook computation. It
does not validate the model against experiment and does not introduce
calibrated physics.

## Condensate Stage

The condensate recovery is closed for the tested deterministic quantities.

Script:

```text
scripts/recover_notebook_condensate_stage.py
```

Outputs:

```text
results/notebook_aligned_recovery/condensate_stage/
```

The stage locks down the notebook defaults in:

```text
configs/notebook_v1_defaults.json
```

Compared quantities include:

- Thomas-Fermi radii;
- chemical potential;
- peak density;
- principal-axis column densities;
- projected Thomas-Fermi profile;
- x-axis projected column-density map;
- central lineouts.

Current deterministic comparison results:

- radii max absolute difference: `0.0`;
- column-density vector max absolute difference: `0.0`;
- projected profile max absolute difference: `0.0`;
- projected column-density map max absolute difference: `0.0`;
- grid shape: `1024 x 1024`;
- peak location: centre pixel `[512, 512]`.

The optional SVG output is tied directly to this recovered quantity and is not
a broad figure-generation workflow.

## Scalar Phase Stage

The scalar phase recovery is closed for the tested deterministic quantities.

Script:

```text
scripts/recover_notebook_phase_stage.py
```

Outputs:

```text
results/notebook_aligned_recovery/phase_stage/
```

Notebook phase references:

- cell 10: `delta_of(...)` and `phi_peak(...)`;
- cell 59: Stage 18.2 phase-map construction;
- cell 67: Step 19.1 scalar phase-map display.

Recovered notebook convention:

```text
delta = 2 * Delta_Hz * 2*pi / Gamma
phi_peak = sigma0 * n_col_peak * delta / (2 * (1 + delta**2))
phase_map = phi_peak * projected_profile
```

Recovered parameters:

- detuning: `1.5e9 Hz`;
- resonant cross section: `7.678673341230136e-14 m^2`;
- natural linewidth: `185353966.5617978 rad s^-1`;
- peak column density: `5.3759624525784675e14 m^-2`;
- dimensionless detuning: `101.69491525423727`;
- peak scalar phase: `0.20294165287929014 rad`.

Current deterministic comparison results:

- dimensionless detuning absolute difference: `0.0`;
- phase-peak absolute difference: `0.0`;
- phase-map max absolute difference: `0.0`;
- phase-map max relative difference: `0.0`;
- phase-map shape: `1024 x 1024`;
- peak location: centre pixel `[512, 512]`.

The optional SVG output, `scalar_phase_stage.svg`, follows the recovered
notebook scalar phase-map quantity only. It does not generate PCI, DGI,
Faraday, camera, or multi-shot figures.

## Regression Tests

Focused regression tests:

```text
tests/regression/test_notebook_condensate_recovery.py
tests/regression/test_notebook_phase_recovery.py
```

These tests check stable numerical outputs and helper/notebook-expression
agreement. They do not test SVG pixel appearance.

## What Remains Uncertain

The recovery still has not locked down end-to-end notebook-aligned recipes for:

- PCI image formation;
- DGI image formation;
- Faraday image formation;
- camera binning/noise workflows;
- multi-shot frame rendering;
- optimisation and SNR operating maps.

Display styling is documented only where it is needed to avoid plotting the
wrong physical quantity. The numerical arrays remain the primary recovery
target.

## Recommended Next Step

The next candidate recovery should be the first phase-dependent imaging path,
most likely one of:

```text
phase_map -> PCI image
phase_map -> DGI image
theta_F map -> Faraday image
```

The next recovery should again begin with canonical notebook defaults,
intermediate numerical comparisons, and only one optional notebook-aligned
figure after the numerical quantity has been matched.
