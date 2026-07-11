# Full Evolving Multishot Accumulated-SNR Analysis

## Purpose

This analysis compares two deliberately separate models:

- an identical-frame clean-loss limit using the same matched-ROI observable;
- an evolving heating- and reabsorption-aware multishot sequence.

The existing peak-pixel clean-loss figure remains unchanged. Its values are
recorded in the model-comparison table for qualitative reference, but they are
not treated as numerically interchangeable with the matched-ROI result.

## Parameters and Provenance

Physical defaults come from `configs/notebook_v1_defaults.json`. The detuning
range and operating point come from `configs/dissertation_plots_v1.json`:

- detuning: 0.75-5.0 GHz on a logarithmic scan including 0.75, 1.0, 1.5, 2.0,
  and 3.0 GHz exactly;
- probe power: 3.5 mW;
- exposure: 40 microseconds;
- imaging axis: x, with the y-z image plane;
- quantum efficiency: 0.4;
- camera bin size: 15 simulation pixels;
- read noise: 7 electrons rms per binned camera pixel;
- condensate-depletion threshold: 30%.

These are Version 1 notebook-aligned parameters and are not experimentally
calibrated.

## Stopping Criterion

The canonical heating model keeps total atom number fixed and reduces the
condensate atom number through recoil heating:

```text
T_next = (T^4 + dE/A_E)^(1/4)
N0 = N_total*[1-(T/Tc)^3]
dE = N_gamma*(1+reabsorption)*E_rec.
```

The threshold variable is

```text
1 - N0(frame)/N0(initial).
```

`N_frames_full` is the integer number of acquired pulses whose post-pulse state
remains at or below 30% condensate depletion. Frame index `i` represents the
pre-pulse state after `i` previous pulses. The first state that crosses the
threshold is excluded. Thus the result has no frame-0/array-length ambiguity.

The analytical clean-loss threshold remains separately named
`N_max_clean` and may be fractional.

## Observable Design

The selected observable is a fixed-spatial-ROI, diagonal-covariance
matched-template SNR.

The ROI is derived from the initial Thomas-Fermi transverse radii with a margin
of two camera pixels. It is fixed across detuning, frame number, and mode. PCI
and DGI therefore use exactly the same spatial support. The image template is
allowed to evolve because the condensate profile and phase evolve.

For each binned camera pixel in the ROI:

```text
signal = (I_atoms - I_reference)*N_ph
variance_shot = I_atoms*N_ph
variance_shot+read = variance_shot + sigma_read^2.
```

The matched-template SNR is

```text
SNR^2 = sum_pixels signal_pixel^2/variance_pixel.
```

This is more defensible than a peak pixel because it collects the distributed
image information and applies the same estimator to both modes. It assumes
independent pixel noise and a known noiseless reference image, matching the
current notebook convention. It is not yet a full calibrated likelihood.

## Framewise Calculation

For every accepted frame, the script:

1. updates temperature from deposited recoil energy;
2. computes the condensate atom number;
3. regenerates the Thomas-Fermi state and radii;
4. recomputes the scalar phase map;
5. propagates the complex scattered field through the configured Fourier
   pupil;
6. forms PCI and finite-OD DGI images;
7. bins to camera pixels;
8. computes shot-only and shot-plus-read matched-template SNR;
9. accumulates `sqrt(sum_i SNR_i^2)`.

The total atom number column is constant in this heating model; condensate atom
number, temperature, phase, image signal, and SNR evolve.

## Frame Budgets

| Detuning (GHz) | Clean-loss continuous `N_max` | Heating continuous stop | Accepted integer frames |
| ---: | ---: | ---: | ---: |
| 0.75 | 7.40 | 3.20 | 3 |
| 1.50 | 29.58 | 13.76 | 13 |
| 3.00 | 118.32 | 56.23 | 56 |

All full-model counts are below the clean-loss continuous threshold. The
integer full-model count has a fitted exponent of 2.121 over 0.75-5 GHz,
compared with 2.000 for clean loss. The larger exponent comes from
detuning-dependent reabsorption plus integer-frame steps; the coefficient is
nevertheless much lower.

## Accumulation Results

The primary full-model quantity is

```text
SNR_total_full = sqrt(sum_i SNR_i^2).
```

The diagnostic identical-frame estimate

```text
SNR_initial*sqrt(N_frames_full)
```

overestimates the explicit result by 8.2%-27.7% across all scanned modes and
noise conditions. The largest errors occur for DGI with read noise because its
signal falls rapidly as the condensate depletes.

At 1.5 GHz:

| Mode | Noise | Full RMS SNR | Identical-frame estimate | Overestimate |
| --- | --- | ---: | ---: | ---: |
| PCI | shot only | 120.33 | 135.22 | 12.4% |
| PCI | shot + read | 117.79 | 132.40 | 12.4% |
| DGI | shot only | 62.82 | 71.09 | 13.2% |
| DGI | shot + read | 21.92 | 26.82 | 22.3% |

## Invariance and Noise Ordering

Over 0.75-3.0 GHz:

- PCI shot-noise accumulated SNR relative range: 10.5%;
- PCI shot-plus-read relative range: 10.2%;
- DGI shot-noise relative range: 4.7%;
- DGI shot-plus-read relative range: 112.7%.

The shot-noise curves therefore retain approximate, not exact, invariance. The
DGI shot-plus-read curve decreases strongly; its high-detuning fitted slope is
`-0.933`, close to the expected read-noise-dominated `-1` accumulated scaling.

For every detuning:

- PCI shot plus read noise is no greater than the PCI shot-noise limit;
- DGI shot plus read noise is no greater than the DGI shot-noise limit;
- the full evolving accumulated SNR is below the same-observable clean-loss
  reference.

## PCI/DGI Comparison

Under the full matched-ROI shot-noise observable, the median PCI/DGI
accumulated-SNR ratio is 1.95. Its relative range across the scan is 19.2%, so
it is not a universal constant.

The median remains close to the legacy peak-pixel factor near two, but the
detuning dependence demonstrates the importance of finite aperture, finite
DGI stop leakage, spatial integration, evolving cloud size, and integer frame
count. The result is estimator-dependent and should not be presented as a
fundamental factor-of-two information advantage.

## Comparison with the Existing Analytical Figure

The same-observable clean-loss matched-ROI curves are an upper bound at every
detuning. Heating and reabsorption reduce the frame count, while evolving
depletion reduces later-frame SNR.

The legacy clean-loss figure uses a peak-pixel analytical observable. The
model-comparison CSV records relative numerical differences, but flags that
those absolute values are not directly comparable. Their qualitative claims
are comparable:

- linear/PCI accumulated SNR remains relatively insensitive to detuning;
- ideal DGI retains approximate photon-limited invariance;
- DGI with read noise falls with detuning;
- clean loss is optimistic.

## Figure Interpretation

Panel (a) shows the identical-frame clean-loss limit using the new matched ROI.
Panel (b) shows the evolving heating sequence with integer stopping and RMS
accumulation. The visible small steps in panel (b) are physical consequences of
integer frame budgets and are not smoothed away.

Caption candidate:

> Accumulated matched-template SNR for PCI and DGI at a nominal 30% condensate-
> depletion budget. Panel (a) uses the identical-frame clean-loss limit with
> the same fixed ROI and image observable; panel (b) accumulates frame-dependent
> SNR as `sqrt(sum_i SNR_i^2)` over the integer frames accepted by the heating-
> and reabsorption-aware sequence. Each frame regenerates the Thomas-Fermi
> state, phase map, Fourier image, and camera-level photon/read variance.
> Results are representative Version 1 outputs and are not experimentally
> calibrated; technical noise, reference-image noise, and density-updated
> reabsorption remain omitted.

## Remaining Limitations

- reabsorption is evaluated from the initial density and held fixed within each
  sequence;
- total atoms do not leave in the heating model;
- reference-image noise is not included;
- pixel covariance and technical probe noise are absent;
- the ROI is fixed geometrically rather than fitted experimentally;
- camera and destruction parameters are uncalibrated.

## Merge Readiness

**A. READY TO MERGE.**

No critical physics discrepancy was found. The full sequence uses the existing
heating/reabsorption equations, explicit integer stopping, evolving
Thomas-Fermi images, matched spatial support, correct RMS accumulation, and
valid same-mode noise ordering. The clean-loss and full-model results remain
clearly separated.
