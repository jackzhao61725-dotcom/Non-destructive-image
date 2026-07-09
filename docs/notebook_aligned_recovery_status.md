# Notebook-Aligned Recovery Status

## Current Recovery Stage

The first canonical notebook-aligned recovery recipe is now:

```text
parameters -> Thomas-Fermi condensate -> projected profile -> column-density map
```

This stage corresponds to the original notebook's condensate construction and
Stage 18.1 density figure. It is intentionally limited to the physical object
used by later phase, imaging, camera, and multi-shot calculations.

## Why Recovery Starts Here

The previous broad figure-gallery attempt diverged because plotting began
before the notebook computational pipeline had been locked down. The cloud
object is upstream of every later quantity:

- scalar phase maps;
- Faraday rotation maps;
- PCI and DGI images;
- camera-level frames;
- multi-shot loss and SNR sequences.

Any notebook-aligned recovery must first prove that the condensate state,
coordinate grid, projected profile, and column-density map match the notebook.

## Canonical Defaults

Notebook Version 1 defaults are collected in:

```text
configs/notebook_v1_defaults.json
```

The config records:

- atom species and mass number;
- atom number;
- scattering length;
- trap frequencies;
- Thomas-Fermi grid size and field of view;
- imaging axis and transverse plane convention;
- unit conversions for `um`, `cm^-3`, and `cm^-2`;
- display metadata for the Stage 18.1 condensate figure;
- a clear note that these are historical notebook defaults, not calibrated
  experimental constants.

## Generated Outputs

The recovery script is:

```text
scripts/recover_notebook_condensate_stage.py
```

Outputs are written to:

```text
results/notebook_aligned_recovery/condensate_stage/
```

Generated files:

- `comparison_report.json`;
- `condensate_summary.json`;
- `central_lineouts.csv`;
- `density_cuts.csv`;
- `metadata.json`;
- `condensate_density_stage.svg`.

The SVG is optional supporting output tied directly to the recovered numerical
quantities. It is not a broad figure-generation workflow and does not generate
phase, PCI, DGI, Faraday, camera, or multi-shot figures.

## What Was Compared

The script compares direct notebook expressions against current helpers:

- `build_thomas_fermi_state(...)`;
- `thomas_fermi_profile_2d(...)`.

Compared quantities:

- chemical potential;
- chemical-potential temperature;
- peak density;
- atom-number consistency check;
- Thomas-Fermi radii;
- principal-axis column densities;
- 2D projected Thomas-Fermi profile;
- x-axis projected column-density map.

Current result:

- radii max absolute difference: `0.0`;
- column-density vector max absolute difference: `0.0`;
- projected profile max absolute difference: `0.0`;
- projected column-density map max absolute difference: `0.0`;
- grid shape: `1024 x 1024`;
- peak location: centre pixel `[512, 512]`.

## Regression Test

The focused regression test is:

```text
tests/regression/test_notebook_condensate_recovery.py
```

It checks stable numerical outputs only. It does not test plotting style or SVG
pixels.

## What Remains Uncertain

This recovery does not validate the notebook physics against experiment. It
only locks down the notebook Version 1 computational prototype for the
condensate stage.

The following stages are not recovered yet:

- scalar phase map;
- Faraday rotation map;
- PCI image;
- DGI image;
- Faraday image;
- camera binning/noise workflows;
- multi-shot frame rendering;
- optimisation and SNR operating maps.

## Recommended Next Step

The next candidate recovery should be the notebook Step 19.1 scalar phase map:

```text
condensate/profile -> phi_peak -> phase_map
```

That recovery should compare:

- `phi_peak(1.5e9, n_col[0])`;
- `phase_map = phi_peak * profile`;
- central lineouts along the notebook `y` and `z` axes;
- display units and axis limits.

No phase or imaging recovery should begin until the condensate-stage outputs
remain stable under tests.
