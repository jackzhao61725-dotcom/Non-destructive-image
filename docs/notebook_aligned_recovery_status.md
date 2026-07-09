# Notebook-Aligned Recovery Status

## Current Recovery Path

The canonical notebook-aligned recovery path now covers:

```text
parameters -> Thomas-Fermi condensate -> projected profile -> column-density map -> scalar phase map -> PCI image
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

Current deterministic comparison results:

- detuning: `1.5e9 Hz`;
- dimensionless detuning: `101.69491525423727`;
- peak scalar phase: `0.20294165287929014 rad`;
- dimensionless detuning absolute difference: `0.0`;
- phase-peak absolute difference: `0.0`;
- phase-map max absolute difference: `0.0`;
- phase-map max relative difference: `0.0`;
- phase-map shape: `1024 x 1024`;
- peak location: centre pixel `[512, 512]`.

## PCI Stage

The PCI recovery is closed for the tested deterministic quantities.

Script:

```text
scripts/recover_notebook_pci_stage.py
```

Outputs:

```text
results/notebook_aligned_recovery/pci_stage/
```

Notebook PCI references:

- cell 16: phase-to-intensity transfer curve and PCI normalisation notes;
- cell 18: Fourier-optics `sim_image(...)` implementation.

Recovered notebook convention:

```text
object_field = np.exp(1j * phase_map)
scattered_field = object_field - 1
propagated_scattered_field = np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)
reference_field = t_p * np.exp(1j * theta)
pci_image_intensity = np.abs(reference_field + propagated_scattered_field) ** 2
```

Recovered parameters:

- imaging axis: `x`, transverse plane `y,z`;
- numerical aperture: `0.08`;
- phase-plate amplitude transmittance: `0.95`;
- phase-plate phase: `pi/2`;
- plate-background intensity: `0.9025`;
- pupil nonzero pixels: `1245`;
- intensity convention: incident-`I0` normalised.

Current deterministic comparison results:

- object-field max absolute difference: `0.0`;
- object-field max relative difference: `0.0`;
- PCI reference-field absolute difference: `0.0`;
- propagated-scattered-field max absolute difference: `1.1102230246251565e-16`;
- propagated-scattered-field max relative difference: `4.117787554919934e-09`;
- PCI image max absolute difference: `0.0`;
- PCI image max relative difference: `0.0`;
- PCI image shape: `1024 x 1024`;
- PCI image peak location: centre pixel `[512, 512]`;
- PCI image centre value: `1.1585677695076864`.

The optional SVG output, `pci_image_stage.svg`, follows the recovered notebook
PCI intensity quantity only. It does not generate DGI, Faraday, camera, or
multi-shot figures.

## Regression Tests

Focused regression tests:

```text
tests/regression/test_notebook_condensate_recovery.py
tests/regression/test_notebook_phase_recovery.py
tests/regression/test_notebook_pci_recovery.py
```

These tests check stable numerical outputs and helper/notebook-expression
agreement. They do not test SVG pixel appearance.

## What Remains Uncertain

The recovery still has not locked down end-to-end notebook-aligned recipes for:

- DGI image formation;
- Faraday image formation;
- camera binning/noise workflows;
- multi-shot frame rendering;
- optimisation and SNR operating maps.

Display styling is documented only where it is needed to avoid plotting the
wrong physical quantity. The numerical arrays remain the primary recovery
target.

## Recommended Next Step

The next candidate recovery should be the DGI imaging path:

```text
phase_map -> DGI image
```

That recovery should reuse the same condensate, phase-map, grid, FFT, and pupil
defaults, then compare only the DGI-specific reference-field and intensity
convention. It should not generate Faraday, camera, or multi-shot figures.
