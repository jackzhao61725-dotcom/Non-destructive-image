# Dissertation Results Generation

## Purpose

This workflow generates dissertation-ready placeholder results for the Version
1 deterministic Faraday optimisation layer.

The outputs are representative and uncalibrated. They are intended to be easy
to regenerate after closed-loop experimental calibration updates the simulator
parameters.

## Files

Configuration:

```text
configs/dissertation_results_v1.json
```

Generator:

```text
scripts/generate_dissertation_results.py
```

Default output directory:

```text
results/faraday_optimisation_v1/
```

## How To Regenerate

From the repository root, run:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
python scripts\generate_dissertation_results.py --config configs\dissertation_results_v1.json
```

The script generates:

- detuning sweep table;
- probe-power sweep table;
- exposure-time sweep table;
- summary JSON;
- metadata JSON;
- simple SVG trade-off figures.

## How To Change Parameters

Adjust parameters in `configs/dissertation_results_v1.json`. Important
configurable entries include:

- `kappa_F`;
- detuning values;
- probe-power values;
- exposure-time values;
- representative density scale;
- cloud-size parameters;
- camera parameters;
- optimisation metric;
- output directory;
- result label.

The script should not be edited for ordinary parameter changes. New
experimental or calibrated results should be produced by creating a new config
file with a new label and output directory.

## Future Calibrated Results

After closed-loop calibration, create a new config such as:

```text
configs/dissertation_results_calibrated_example.json
```

That calibrated config should record the updated density scale, optical-depth
scale, camera parameters, and any experimentally fitted Faraday calibration
parameters.

Dissertation figures generated from the current Version 1 config should be
treated as placeholders. They may be replaced after absorption / RAI calibration
and future `kappa_F` fitting.

## Metadata Policy

The generated `metadata.json` states:

- the result label;
- that the outputs are Version 1 representative / uncalibrated results;
- the configured `kappa_F`;
- that no experimental RAI / absorption calibration has yet been applied;
- that outputs should be regenerated after closed-loop calibration;
- the git commit hash used when generating results;
- the config file used;
- the helper functions used.

## Current Limitations

The current result outputs are:

- deterministic only;
- not noise averaged;
- not fitted to RAI or absorption data;
- still based on the Version 1 phenomenological Faraday convention;
- not based on experimentally fitted `kappa_F`;
- not based on a microscopic Faraday model;
- not a substitute for final calibrated dissertation figures.

They are suitable as Version 1 representative outputs and as a reproducible
starting point for later calibrated result generation.
