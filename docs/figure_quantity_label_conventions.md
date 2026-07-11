# Figure Quantity Label Conventions

## Purpose

Dissertation figures must label the physical quantity being plotted, not the
implementation variable that produced it. If a precise subscript is not defined
in the dissertation text, use a generic label rather than inventing new
notation.

This document applies to notebook-aligned recovery figures, representative V1
plots, numerical audit plots, and future calibrated figures.

## General Rules

- Figure labels must match the dissertation's physical definitions.
- Panel titles should name the plotted quantity or comparison, not the workflow
  stage that generated it.
- Captions and metadata should carry parameter details, caveats, and exact
  formulas.
- Do not call a signal, phase, rotation, density, count, or fitted residual an
  intensity.
- If plotted values are normalised, metadata must state what was normalised and
  by what reference.
- If absolute values are saved in CSV/JSON alongside a normalised plot, metadata
  must say so.

## Intensity

Use:

- `I` for generic intensity;
- `I/I0` for intensity normalised to incident probe intensity;
- `camera counts` only when the plotted quantity is actually electron or photon
  counts;
- `normalised intensity` when the exact incident-intensity reference is less
  explicit.

Avoid:

- mode-specific intensity subscripts unless they are defined in the text;
- using `intensity` for dimensionless difference signals such as `S`;
- using camera-count language for values already divided by the photon scale.

## Dual-Port Faraday

Use:

- `I_H` and `I_V` for the two analysed linear-polarisation port intensities;
- `S` for the normalised difference signal.

Definition:

```text
S = (I_H - I_V) / (I_H + I_V)
```

`S` is a signal, not an intensity. Colourbars, lineouts, captions, summaries,
and metadata must not label `S` as `I/I0`, intensity, or camera counts.

For notebook-aligned recovery scripts that still use historical notebook port
names such as `u` and `v`, metadata should record the mapping to dissertation
notation. For the current Faraday camera panel:

```text
I_H corresponds to the notebook v port.
I_V corresponds to the notebook u port.
S = (I_H - I_V) / (I_H + I_V)
  = (I_v - I_u) / (I_v + I_u)
```

This preserves the notebook numerical convention while presenting the
dissertation-facing notation cleanly.

## Dark-Field Faraday

Use:

- `I/I0` when the plotted quantity is normalised transmitted intensity;
- `normalised intensity` when the exact incident-intensity convention is not
  explicit.

Do not introduce special dark-field subscripts in figures unless the
dissertation text defines them. It is acceptable for metadata to mention the
notebook variable `I_dark`, but the visible figure label should remain a
standard intensity label unless a special symbol is defined.

## PCI and DGI

Use:

- `I/I0` for normalised image intensity when the notebook/helper output is an
  incident-intensity-normalised image;
- `normalised intensity` when the exact reference is less explicit.

Avoid unclear mode-specific subscripts in figure labels unless they have been
defined in the surrounding dissertation text.

## Phase and Faraday Rotation

Use:

- `phi` for scalar optical phase, with units of radians;
- `theta_F` for Faraday rotation angle, with units of radians.

Do not call `phi` or `theta_F` an intensity or image signal. If a plotted
quantity is derived from them, label the derived quantity instead, for example
`S` for the dual-port difference.

## Column Density

Use:

- `n_col(y,z)` for the full 2D column-density distribution in the `y,z` plane;
- `n_col(x,z)` for the full 2D column-density distribution in the `x,z` plane;
- `n_col(x,y)` for the full 2D column-density distribution in the `x,y` plane.

Reserve:

- `tilde n_x`, `tilde n_y`, and `tilde n_z` for peak column-density scalar
  values, where the subscript denotes the line-of-sight integration axis.

Do not use `tilde n_i` for a full 2D distribution.

## Scattering and Destructiveness

Use:

- `N_gamma` for scattered photons per atom per exposure;
- `R_sc` for scattering rate.

Do not mix them unless the figure or caption explicitly states:

```text
N_gamma = R_sc * tau
```

where `tau` is the exposure time.

## Normalisation Metadata

Every generated figure with normalised plotted values should record:

- plotted quantity;
- symbol used in the figure;
- reference value or reference image;
- whether plotted values are absolute, normalised, binned, noisy, or counts;
- whether absolute values are saved elsewhere;
- calibration status.

For notebook-aligned V1 figures, metadata should also state that the output is
uncalibrated and that no experimental RAI/absorption calibration has been
applied.

## Current Faraday Camera Panel Labels

The Faraday camera panel uses:

- dark-field image colourbar: `I/I0`;
- dark-field lineout y-axis: `I/I0`;
- dual-port signal map colourbar: `S`;
- dual-port lineout y-axis: `S`;
- port intensities in metadata: `I_H/I0` and `I_V/I0`;
- signal definition in metadata: `S = (I_H - I_V)/(I_H + I_V)`.

This makes clear that `S` is not an intensity.

