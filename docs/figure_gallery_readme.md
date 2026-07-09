# Dissertation Figure Gallery V1

## Purpose

The Version 1 figure gallery provides dissertation-facing representative
figures generated from the current modular helper package.

These figures are placeholders for writing and presentation. They are not final
experimental results, not calibrated predictions, and not copied directly from
the exploratory notebook.

The original notebook is used only as historical context for candidate
calculations and visualisations. The generated figures use the current modular
helpers where possible.

## Regeneration

From the repository root, run:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
python scripts\generate_dissertation_figures.py --config configs\dissertation_figures_v1.json
```

Outputs are written to:

```text
results/dissertation_figures_v1/
```

## Configuration

Adjust figure parameters through:

```text
configs/dissertation_figures_v1.json
```

The config contains:

- output directory;
- `kappa_F`;
- representative density and cloud-size parameters;
- detuning, probe power, and exposure time;
- one-dimensional optimisation sweep values;
- camera parameters;
- fixed RNG seed for the noisy camera figure;
- multi-shot sequence parameters;
- synthetic absorption / RAI calibration-readiness parameters;
- notes that no experimental absorption / RAI calibration has been applied.

Create a new config and output directory for future calibrated figures rather
than overwriting the Version 1 representative gallery.

## Generated Figures

The gallery is intentionally kept in one directory rather than split into
main-text, appendix, and backup layers. Use the figure index to decide which
figures are most relevant to a dissertation section.

- `closed_loop_architecture.svg`
- `tf_density_model.svg`
- `tf_profile_lineouts.svg`
- `light_atom_tradeoff.svg`
- `phase_and_faraday_maps.svg`
- `imaging_mode_comparison.svg`
- `imaging_mode_lineouts.svg`
- `camera_realism.svg`
- `multishot_evolution.svg`
- `optimisation_sweep_summary.svg`
- `absorption_calibration_observables.svg`

Supporting data and metadata:

- `figure_index.json`
- `tf_density_model_data.csv`
- `tf_profile_lineouts_data.csv`
- `light_atom_tradeoff_data.csv`
- `phase_and_faraday_maps_metadata.json`
- `imaging_mode_comparison_metadata.json`
- `imaging_mode_lineouts_data.csv`
- `camera_realism_metadata.json`
- `multishot_evolution_data.csv`
- `optimisation_sweep_summary_data.csv`
- `optimisation_sweep_summary_metadata.json`
- `absorption_calibration_observables_data.csv`
- `absorption_calibration_observables_metadata.json`
- `summary.json`
- `metadata.json`

`figure_index.json` records, for every SVG figure, what it shows, which helpers
or script logic were used, whether it is deterministic, and whether it depends
on placeholder Version 1 parameters.

## Current Limitations

- Figures are Version 1 representative / uncalibrated outputs.
- `kappa_F = 1.0` remains a placeholder where Faraday figures are involved.
- No experimental RAI / absorption calibration has been applied.
- `kappa_F` has not been experimentally fitted.
- No microscopic Faraday model is introduced.
- Camera realism uses one fixed-seed noisy frame, not stochastic averaging.
- Multi-shot evolution is deterministic bookkeeping, not full noisy frame
  rendering.
- The optimisation figure shows separate one-dimensional sweeps, not a full
  multi-parameter optimisation.
- The absorption / RAI calibration figure is synthetic and demonstrates
  observable extraction only; it is not fitted to laboratory data.

## Overclaiming Risks

Use these figures to explain the simulator and representative trade-offs. Do
not present them as final experimental operating-point predictions.

After closed-loop calibration, regenerate figures using calibrated density,
camera, optical-depth, and Faraday parameters in a new config file.
