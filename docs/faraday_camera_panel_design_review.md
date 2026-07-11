# Faraday Camera-Level Reference Figure Design Review

## Scope

This document reviews whether the next notebook-aligned dissertation figure
should be a Faraday camera-level reference panel. It is a design review only:
no figure, SVG, PNG, simulator change, helper API change, or physics change is
introduced here.

The intended figure family is:

```text
ideal Faraday fields -> camera binning/noise -> dark-field and dual-port
camera-level reference panel
```

The correct readiness classification is:

```text
READY AFTER SMALL ORCHESTRATION
```

The figure is scientifically useful, but it should not be generated until a
small notebook-aligned Faraday camera orchestration script exists and has
written comparison metadata.

## Notebook Counterpart

The primary notebook counterpart is `notebook_sections/06_faraday.py` cell 51,
lines 102-170 in the exported file.

That cell defines:

- `sim_faraday_fields(axis, theta_F_val)`;
- `faraday_maps(Delta_Hz, axis=0)`;
- a reference operating point with `Delta = 1.5e9 Hz`, `axis = 0`;
- camera-level dark-field and dual-port outputs using `to_camera(...)`;
- a 2 x 2 figure saved in the notebook as `fig_faraday_reference.png`.

What cell 51 shows:

- top left: noisy camera-level dark-field Faraday image,
  `cam_dark = to_camera(fm["I_dark"], 5.0)[0]`;
- top right: noisy camera-level dual-port normalised difference,
  `S_map = (cam_v - cam_u) / (cam_v + cam_u)`;
- bottom left: central dark-field lineout comparing noiseless binned image
  `ideal_dark` against noisy camera pixels `cam_dark`;
- bottom right: central dual-port lineout comparing noiseless binned
  `S_ideal` against noisy camera-level `S_map`.

It includes both dark-field and dual-port Faraday. It is not only an ideal
field figure and not only a port-intensity figure.

The camera/noise recipe comes from `to_camera(...)` in
`notebook_sections/07_camera.py` cell 20:

```text
bin high-resolution I/I0 image by 15 x 15 averaging
Nd = N_phot_pix(P_mW, tau_s=None, QE=None)
counts = rng.poisson(clip(binned, 0, None) * Nd) + rng.normal(0, read_e, shape)
camera_image = counts / Nd
```

The notebook RNG source is `notebook_sections/00_imports.py` cell 0:

```text
rng = np.random.default_rng(7)
```

Important RNG caveat:

- cell 51 depends on notebook-global RNG state because it calls `to_camera`
  three times for `I_dark`, `I_u`, and `I_v`;
- exact pixel replay is only meaningful if the future recovery script either
  reproduces the relevant call order explicitly or defines a clean explicit
  seed replay and reports it as such;
- a future figure should not claim exact arbitrary interactive notebook RNG
  reproduction.

Related notebook cells:

- cell 85: ideal Faraday rotation / `|E_y|^2` explanatory panel;
- cell 87: dark-field Faraday recorded panel with noiseless, noisy, and lineout
  views;
- cell 89: ideal dual-port `I_u`, `I_v`, and `S` panel;
- cell 91: dual-port noisy port frames plus flicker-cancellation demonstration.

Those related cells are useful context, but the clean next recovery target
should be cell 51. Cell 91 should remain a later, separate figure because it
adds flicker modelling and a different narrative point.

## Current Recovery Support

Already available:

- canonical notebook defaults in `configs/notebook_v1_defaults.json`;
- recovered condensate stage under
  `results/notebook_aligned_recovery/condensate_stage/`;
- recovered scalar phase stage under
  `results/notebook_aligned_recovery/phase_stage/`;
- recovered ideal Faraday stage under
  `results/notebook_aligned_recovery/faraday_stage/`;
- deterministic camera recovery under
  `results/notebook_aligned_recovery/camera_stage/`;
- stochastic camera recovery under
  `results/notebook_aligned_recovery/noisy_camera_stage/`;
- linear approximation audit under `results/linear_approximation_audit/`;
- unit and parameter inventories under `results/notebook_aligned_recovery/`.

Current helpers that can be reused:

- `simulate_faraday_image(...)` for `theta_f_map`, sigma fields, `Ex`, `Ey`,
  dark-field intensity, dual-port `I_u`, `I_v`, and `S`;
- `simulate_camera_image(...)` for deterministic binning and normalisation;
- `simulate_noisy_camera_image(...)` for explicit-seed stochastic camera
  replay;
- `bin_to_camera_pixels(...)`, `add_camera_noise(...)`, and
  `normalize_camera_counts(...)` for direct composition when separate port
  handling is clearer.

Missing small orchestration:

- a Faraday-specific camera recipe that applies the camera/noise model to
  `I_dark`, `I_u`, and `I_v` using the cell 51 operating point;
- construction of `S_map` from noisy port frames and `S_ideal` from
  deterministic binned port frames;
- comparison metadata for binned/noisy dark-field and dual-port arrays;
- lineout extraction using notebook camera-axis convention.

No new public `src/` helper is required for this design. A script/config
recovery step is enough.

## Parameter Provenance

No figure should be generated without recording these parameters.

| Parameter | Value / convention | Provenance | Status |
| --- | --- | --- | --- |
| species | `166Er` | `configs/notebook_v1_defaults.json` / notebook parameters | inherited |
| transition wavelength | `401e-9 m` | `configs/notebook_v1_defaults.json` | inherited |
| natural linewidth | `2*pi*29.5e6 rad/s` | `configs/notebook_v1_defaults.json` | inherited |
| resonant cross section | `3*lambda^2/(2*pi)` | derived from notebook constants | derived |
| atom number | `25000` | `configs/notebook_v1_defaults.json` | inherited |
| scattering length | `72 a0` | `configs/notebook_v1_defaults.json` | inherited |
| trap frequencies | `[293, 14, 233] Hz` | `configs/notebook_v1_defaults.json` | inherited |
| grid | `1024 x 1024`, `FOV = 100 um` | `configs/notebook_v1_defaults.json` | inherited |
| imaging axis | `axis = 0`, displayed as `y,z` plane | notebook cell 51 / config imaging geometry | inherited |
| detuning | `1.5e9 Hz` | notebook cell 51 and phase/Faraday recovery | stage-specific |
| `kappa_F` | `1.0` | notebook cell 49 / config Faraday recovery | inherited placeholder |
| numerical aperture | `0.08` | notebook imaging setup / recovered Faraday metadata | inherited |
| probe power | `5.0 mW` | notebook cell 51 `to_camera(..., 5.0)` | stage-specific |
| exposure time | `100 us` default | notebook camera helper `N_phot_pix(..., tau_s=None)` | inherited camera default |
| camera QE | `0.40` | notebook parameters / config camera recovery | inherited |
| camera binning | `15 x 15`, output `68 x 68` | notebook camera helper / recovered camera stage | inherited |
| read noise | `7 e- rms / pixel` | notebook parameters / config camera recovery | inherited |
| RNG seed | `7` | notebook import cell global RNG | inherited, but call-order sensitive |
| dark-field display scale | `vmin=-0.005`, `vmax=0.05`, `inferno` | notebook cell 51 | display-only |
| dual-port display scale | `vmin=-0.5`, `vmax=0.5`, `RdBu_r` | notebook cell 51 | display-only |
| lineout axis | camera central row against `ycam` in `um` | notebook cell 51 | display-only |

Ambiguities to resolve in the future generation script:

- whether to replay the exact cell 51 RNG call order from a clean seed or to
  label the figure as an explicit-seed replay;
- whether to preserve notebook's cell 51 colour limits exactly or use the
  standardised dissertation label conventions while keeping the same physical
  display ranges;
- whether the future config should add a dedicated
  `faraday_camera_panel_recovery` block, because the existing
  `camera_recovery.probe_power_mw = 2.0` belongs to the PCI camera stage, while
  cell 51 uses `5.0 mW`.

## Scientific Validity

The planned figure should be described as finite-phase / finite-rotation
numerical propagation followed by notebook camera binning/noise. It should not
be described as a small-angle-only plot.

Relevant linear-approximation audit results:

- canonical peak scalar phase: `0.202941652879 rad`;
- canonical peak Faraday rotation with `kappa_F = 1.0`:
  `0.202941652879 rad`;
- dark-field `sin^2(theta_F)` versus `theta_F^2` differs by about `1.384%`
  at the peak;
- dual-port `sin(2 theta_F)` versus `2 theta_F` differs by about `2.799%`
  at the peak.

Correct interpretation:

- the simulation uses exact `exp(+/- i theta_F)` propagation, exact `Ex/Ey`
  recombination, exact `I_dark = |Ey|^2`, and exact dual-port ratio
  `(I_v - I_u)/(I_v + I_u)`;
- small-angle formulae are useful only for scaling explanation in the caption
  or text;
- `kappa_F = 1.0` is a Version 1 phenomenological placeholder;
- no microscopic Faraday model is introduced;
- no experimental RAI/absorption calibration has been applied;
- the figure is not a final operating-point prediction.

Overclaim risks:

- calling the panel calibrated experimental Faraday imaging;
- implying `kappa_F = 1.0` has been fitted to Er atomic-structure data;
- presenting the 5 mW / 100 us camera frame as an optimised operating point;
- saying the response is strictly linear rather than exact finite rotation with
  a weak-rotation interpretation;
- treating the noisy frame as a unique notebook-global RNG replay unless the
  call order is explicitly reproduced.

## Usefulness Decision

Decision:

```text
Generate this figure after a small orchestration script is approved.
```

Classification:

```text
READY AFTER SMALL ORCHESTRATION
```

Reason:

- it fills a documented recovery gap: ideal Faraday outputs have been
  recovered, and generic camera/noise stages have been recovered, but the
  Faraday-specific camera panel has not;
- it is directly tied to notebook cell 51, not an invented broad gallery;
- it is useful for dissertation narrative because it shows the practical
  difference between dark-field and dual-port Faraday at the camera level;
- it bridges the theory/ideal-image stage and the measurement-realism stage;
- it should be main-text candidate material if visually clear, otherwise an
  appendix figure.

It should not wait for experimental calibration, provided it is labelled as
notebook-aligned Version 1 / uncalibrated. Calibration would later replace or
supplement the parameter values, not block recovery of the historical prototype
figure.

## Recommended Figure Design

Recommended layout:

```text
2 x 2 panel, matching notebook cell 51
```

Panels:

1. noisy camera-level dark-field Faraday image;
2. noisy camera-level dual-port signal map `S`;
3. dark-field central lineout: noiseless binned versus noisy camera frame;
4. dual-port central lineout: noiseless binned `S` versus noisy camera-level
   `S`.

Quantities and units:

- image axes: `y (um)` and `z (um)`;
- dark-field image and lineout: `I_dark/I0`;
- dual-port map and lineout: `S = (I_v - I_u)/(I_v + I_u)`;
- use camera-pixel binned arrays for plotted camera-level maps;
- keep ideal high-resolution Faraday fields in metadata/comparison, not as
  extra panels.

Normalisation:

- dark-field and port outputs are already incident-`I0` normalised;
- dual-port `S` is dimensionless;
- no additional normalisation should be introduced for the figure.

Dark-field and dual-port together:

- show both in one figure because cell 51 is explicitly a comparison panel;
- do not include dual-port flicker cancellation in this figure, because that
  belongs to notebook cell 91 and would overcrowd the reference panel.

Figure versus caption:

- figure should show image shape and lineout contrast;
- caption should carry `Delta = 1.5 GHz`, `P = 5 mW`, `tau = 100 us`,
  `QE = 0.40`, `read_e = 7 e-`, `kappa_F = 1.0 placeholder`, and the
  uncalibrated Version 1 caveat;
- caption should state that exact finite-rotation propagation was used and
  small-angle formulae are only interpretive.

## Proposed Follow-Up Task

Do not execute this task until approved.

```text
Task: Recover Notebook-Aligned Faraday Camera-Level Reference Panel

Work from:
work/faraday-camera-panel-design or a new approved feature branch.

Allowed changes:
- configs/notebook_v1_defaults.json only if adding a dedicated display/config
  block for the Faraday camera panel is necessary;
- scripts/recover_notebook_faraday_camera_panel.py;
- results/notebook_aligned_recovery/faraday_camera_panel/;
- docs/notebook_aligned_recovery_status.md;
- tests/regression/test_notebook_faraday_camera_panel.py if stable.

Do not modify:
- src/;
- helper APIs;
- simulator physics;
- notebook_sections/;
- original notebook;
- existing baselines;
- existing recovery data except new faraday_camera_panel outputs;
- dissertation plot data;
- deliverables zip.

Requirements:
1. Reproduce notebook_sections/06_faraday.py cell 51.
2. Use recovered condensate, phase, and ideal Faraday logic.
3. Apply notebook camera binning/noise to I_dark, I_u, and I_v.
4. Construct S_map from noisy port frames and S_ideal from deterministic
   binned port images.
5. Save comparison_report.json, faraday_camera_panel_summary.json,
   central_lineouts.csv, frame_statistics.csv, and metadata.json.
6. Generate only faraday_camera_reference_panel.svg.
7. State RNG policy explicitly.
8. Run pytest, notebook section validation, and the new recovery script.
```

