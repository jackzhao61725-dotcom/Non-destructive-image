# Accumulated-SNR Invariance Physics Audit

For the complete physics-first and code-traceability review, including
clean-loss versus heating/reabsorption benchmarks and integer-frame semantics,
see `docs/accumulated_snr_full_physics_audit.md`.

## Scope

This audit checks the accumulated-SNR figure against the current Version 1
notebook-equivalent implementation. The figure is intended to explain scaling,
not to replace the full Fourier-imaging, camera, or multishot calculations.

## Conclusion

The central scaling argument is physically consistent with the current
notebook in its far-detuned, small-phase analytical regime:

- scalar phase: `|phi| ~ 1/|Delta|`;
- scattered photons per atom: `N_gamma ~ 1/Delta^2`;
- fixed clean-loss budget: `N_max ~ 1/N_gamma ~ Delta^2`;
- accumulated SNR: `SNR_total = SNR_shot * sqrt(N_max)`.

The earlier label `realistic` was too broad. Those curves included photon and
read-noise variance, but did not run the full spatial Fourier/pupil, camera
binning, and sequence-depletion pipeline. They are now labelled `shot + read
noise`.

## Mode Scaling

For PCI, the notebook small-phase phase-plate contrast is
`Delta I = 2 t_p phi`. Photon noise is set mainly by the bright phase-plate
background. Both the photon-shot-noise-only and shot-plus-read-noise per-frame
SNR therefore scale approximately as `1/|Delta|`; multiplying by
`sqrt(N_max) ~ |Delta|` gives an approximately flat accumulated SNR.

For ideal DGI in the opaque-stop limit, the detected intensity is
`4 sin^2(phi/2)`, which approaches `phi^2`. In the photon-shot-noise limit its
per-frame SNR is proportional to `sqrt(I_DGI) ~ |phi| ~ 1/|Delta|`, so the
accumulated SNR is also approximately flat. When read noise dominates, the
per-frame SNR instead follows `I_DGI/read_noise ~ 1/Delta^2`; the accumulated
SNR then falls approximately as `1/|Delta|`.

## Cross-Mode Prefactor

The PCI and DGI shot-noise-limit curves do not coincide. Across the configured
scan, the PCI/DGI accumulated-SNR ratio is about 2.00. This is expected from the
current observable definitions:

- PCI uses the carrier as a homodyne reference and gives approximately
  `2|phi| sqrt(N_ph)`;
- DGI directly counts dark-field photons and gives approximately
  `|phi| sqrt(N_ph)`.

The factor is therefore a mode-dependent measurement prefactor, not a hidden
post-processing normalisation. It should not be interpreted as a universal
factor-of-two experimental advantage because finite aperture, spatial
integration, stop leakage, and fitted backgrounds can change the quantitative
comparison.

## Destruction Model Boundary

The plotted `N_max` uses the notebook clean-loss expression
`-log(1-f)/(eta N_gamma)` at `f = 0.30`. The notebook later identifies this as
an optimistic bound and treats heating plus reabsorption as the more realistic
deep-trap model. Therefore the figure demonstrates the clean analytical
invariance argument; it is not a quantitative prediction of the final usable
frame count under the later destruction model.

The mode independence of `N_max` is still correct: PCI and DGI modify the light
after the atom-light interaction, so the scattering budget upstream is common
to both modes for the same probe power, duration, and detuning.

## Numerical Checks

The generated summary records the executable checks. At the current settings:

- maximum phase is below the notebook `phi < 0.5` linear-regime threshold;
- `N_max` has a fitted detuning exponent close to `+2`;
- PCI accumulated SNR is flat to substantially better than 1%;
- DGI shot-noise-limit accumulated SNR is flat to within 1%;
- DGI with shot plus read noise decreases with detuning;
- shot-plus-read SNR never exceeds the corresponding shot-noise limit;
- the PCI/DGI ideal ratio is close to the expected factor of two.

## Not Included

The figure does not include:

- finite-NA Fourier propagation or axis-dependent blur;
- camera binning or resolution-element integration;
- the notebook DGI stop leakage at `OD = 4`;
- frame-by-frame condensate depletion;
- heating plus reabsorption in the destruction budget;
- experimental camera or absorption/RAI calibration.

These omissions are acceptable for a scaling figure only when stated in the
caption and metadata. A future quantitative operating-point figure should use
the full recovered imaging, camera, and multishot pipeline.
