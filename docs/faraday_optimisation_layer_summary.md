# Version 1 Faraday Optimisation Layer Summary

## Status

The repository now includes a completed deterministic single-variable Faraday
optimisation layer above the Version 1 migrated simulator core.

This layer is documentation- and regression-test-supported, but it remains
conservative: it does not redesign the Faraday model, change helper APIs, alter
notebook logic, or recalibrate `kappa_F`.

## Completed Helpers

The public optimisation helpers are:

- `evaluate_faraday_operating_point(...)`
- `sweep_faraday_detuning(...)`
- `sweep_faraday_intensity(...)`
- `sweep_faraday_exposure_time(...)`
- `summarise_faraday_sweep(...)`

`evaluate_faraday_operating_point(...)` evaluates one deterministic Faraday
operating point and reports signal, scattering, destructiveness, and
information-efficiency metrics.

The three sweep helpers evaluate one operating variable at a time while holding
the other operating parameters fixed.

`summarise_faraday_sweep(...)` converts a sweep result into a small report-ready
dictionary containing the selected metric, best index, best parameter value,
best metric value, and number of evaluated points.

## Optimisation Philosophy

The optimisation target is not maximum image brightness alone.

For continuous non-destructive imaging, a bright image can be a poor operating
point if it causes excessive scattering, heating, reabsorption, or condensate
loss. The relevant question is how much useful information is gained for a
given destructive cost.

The current deterministic metrics therefore include quantities such as:

- `signal_per_scattered_photon`
- `information_per_scattered_photon`
- `signal_to_destruction`
- `estimated_per_frame_snr`

These metrics connect image signal to the information-versus-destruction
research question.

## One-Dimensional Sweeps

### Detuning Sweep

`sweep_faraday_detuning(...)` evaluates the operating-point objective across a
list of probe detunings.

This is useful because the dispersive Faraday signal and the scattering cost
scale differently with detuning. The best information-efficiency point may not
be the brightest image.

### Probe Intensity Sweep

`sweep_faraday_intensity(...)` evaluates the objective across a list of
`probe_power_mw` values.

This provides a deterministic way to compare signal efficiency and destructive
cost as the probe intensity-like control parameter changes, without introducing
a multi-parameter search.

### Exposure-Time Sweep

`sweep_faraday_exposure_time(...)` evaluates the objective across a list of
`pulse_duration_s` values.

This supports deterministic exposure-time scaling studies while keeping the
current Faraday model and lower-layer helper behaviour unchanged.

## Summary Helper

`summarise_faraday_sweep(...)` accepts a deterministic sweep result and a metric
name. It selects the best parameter value according to that metric and reports:

- metric name;
- parameter name;
- best index;
- best parameter value;
- best metric value;
- number of evaluated points;
- small min/max diagnostics for available quantities.

It does not plot, write files, use pandas, run a new sweep, or perform
stochastic averaging.

## Limitations

The current Faraday optimisation layer does not yet include:

- two-dimensional or three-dimensional sweeps;
- plotting or operating-map generation;
- stochastic noise averaging;
- experimental RAI calibration;
- `kappa_F` recalibration;
- microscopic Faraday physics;
- automated optimisation algorithms.

The current Faraday convention remains:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

`kappa_F` remains fixed and phenomenological in this Version 1 optimisation
layer.

## MSc Relevance

This layer provides the first optimisation-ready interface above the migrated
simulator core. It allows the MSc project to move from image simulation toward
quantitative operating-point comparison while preserving the notebook-equivalent
physics.

The layer supports future optimisation over detuning, probe intensity, and
exposure time, and it frames those choices in terms of information gained versus
destructiveness. This makes it directly relevant to continuous
non-destructive imaging of an ultracold `166Er` Bose-Einstein condensate.

## Next Directions

Future work should remain staged:

- add report-ready tables or plotting as a presentation layer;
- design two-dimensional sweeps before implementing them;
- add stochastic noise averaging only with explicit RNG policy;
- add RAI-based calibration in a separate calibration layer;
- revisit `kappa_F` only in a calibration-specific milestone.
