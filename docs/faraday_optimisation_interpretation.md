# Faraday Optimisation Interpretation Notes

## Optimisation Goal

The Faraday optimisation objective is not to maximise image brightness alone.
The project goal is to maximise useful information gained from the condensate
while minimising destructiveness.

In this context, a brighter image is only useful if the extra signal is worth
the additional scattering, heating, reabsorption, atom loss, or depletion of the
available photon budget. A dimmer image may be preferred if it gives more
useful information per scattered photon.

## Role of Detuning

Probe detuning controls both dispersive signal and destructiveness.

For the current Version 1 model:

- the Faraday / dispersive signal decreases as the probe is moved farther from
  resonance;
- scattering and destructiveness decrease faster with detuning in the
  far-detuned limit;
- therefore, a useful information-versus-destruction metric may favour a
  detuning away from the brightest image.

This is why the best detuning from the current sweep should not be interpreted
as the detuning that gives the largest raw Faraday signal. It is the detuning
that scores best under the chosen objective metric.

## Supported Metrics

`evaluate_faraday_operating_point(...)` currently returns deterministic scalar
quantities for one operating point:

- `faraday_signal_rad`
- `faraday_signal_scale`
- `scattered_photons_per_atom`
- `reabsorption_fraction`
- `destructiveness_metric`
- `estimated_per_frame_snr`
- `signal_per_scattered_photon`
- `information_per_scattered_photon`
- `signal_to_destruction`

`sweep_faraday_detuning(...)` evaluates the same quantities for a list of
detunings and stacks them into arrays. It also reports:

- `objective_key`
- `best_index`
- `best_detuning_hz`
- `best_objective_value`

The default sweep objective is `signal_per_scattered_photon`.

## Interpreting Best Detuning

When `objective_key="signal_per_scattered_photon"`, the reported
`best_detuning_hz` identifies the detuning with the largest Faraday signal scale
per scattered photon. This means:

- it is an efficiency optimum under the current deterministic model;
- it is not necessarily the brightest image;
- it does not yet include stochastic camera-noise averaging;
- it should be treated as a candidate operating point, not a final experimental
  recommendation.

Other supported metrics can be selected if the scientific question changes. For
example, `signal_to_destruction` includes the current destructiveness estimate,
while `estimated_per_frame_snr` emphasises expected per-frame image quality when
a photon count scale is supplied.

## Limitations

The current optimisation helpers are deliberately small. They do not yet
include:

- intensity sweep;
- exposure-time sweep;
- stochastic noise averaging;
- full camera-noise objective;
- multi-parameter optimisation;
- plotting or operating-map generation;
- experimental RAI calibration;
- recalibration of the phenomenological Faraday factor.

The current Faraday convention remains:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

`kappa_F` is phenomenological and fixed during the current migration and early
optimisation work. It should only be revisited in a later calibration-specific
milestone.

## Recommended Next Milestone

After confirming that the detuning sweep metrics behave sensibly, the next
implementation milestone should add one small additional sweep dimension:

- either a probe-intensity sweep; or
- an exposure-time scaling sweep.

That milestone should remain deterministic and should still avoid plotting,
stochastic averaging, or full multi-parameter optimisation until the objective
definitions are stable.
