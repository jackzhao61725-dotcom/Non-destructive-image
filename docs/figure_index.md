# Figure and data index

This index lists the maintained dissertation-facing outputs. Each generator
writes numerical data and metadata alongside the rendered figure.

## Active figures

| Output | Generator | Config | Quantity and scope |
| --- | --- | --- | --- |
| `results/dissertation_plots_v1/accumulated_snr_invariance/figure_3_2.{png,pdf}` | `scripts/generate_accumulated_snr_invariance_plot.py` | `configs/dissertation_plots_v1.json` | analytical peak-pixel accumulated SNR under the continuous clean-loss budget |
| `results/dissertation_plots_v1/figure_4_2/figure_4_2.{png,pdf,svg}` | `scripts/generate_figure_4_2.py` | `configs/figure_4_2.json` | noiseless finite-aperture PCI, DGI, dark-field Faraday and dual-port outputs |
| `results/dissertation_plots_v1/figure_5_1/figure_5_1.{png,pdf,svg}` | `scripts/generate_figure_5_1.py` | `configs/figure_5_1.json` | single-frame `SNR_3x3` fluence scan and representative noisy Faraday frames |
| `results/dissertation_plots_v1/figure_5_2/figure_5_2_dual_port_heatmap.{png,pdf,svg}` | `scripts/generate_figure_5_2.py` | `configs/figure_5_2.json` | strict usable dual-port frame count over detuning and fluence |
| `results/dissertation_plots_v1/figure_5_2/figure_5_2_operating_band.{png,pdf,svg}` | `scripts/generate_figure_5_2.py` | `configs/figure_5_2.json` | fixed-1.5-GHz operating-band comparison |
| `results/dissertation_plots_v1/figure_5_4/figure_5_4.{png,pdf,svg}` | `scripts/generate_figure_5_4_snr_panel.py` | `configs/figure_5_4.json` | 25-frame dual-port sequences at 30, 50, 90 and 150 mW us |
| `results/dissertation_plots_v1/full_multishot_accumulated_snr/full_multishot_accumulated_snr.{png,pdf,svg}` | `scripts/generate_full_multishot_accumulated_snr.py` | `configs/dissertation_plots_v1.json` | four-mode evolving matched-template accumulated SNR |

The `figure_5_2` directory also contains a combined overview figure and the
complete scan CSV/JSON. Which rendered panel receives the final dissertation
figure number is a writing decision; the numerical source remains one scan.

## Supporting physics plot

| Output | Generator | Config | Scope |
| --- | --- | --- | --- |
| `results/dissertation_plots_v1/detuning_tradeoff/detuning_tradeoff.svg` | `scripts/generate_detuning_tradeoff_plot.py` | `configs/dissertation_plots_v1.json` | normalised dispersive-signal and scattering trends |

## Historical notebook-aligned figures

The directories below are retained as regression and interpretation evidence,
not as current detector predictions:

- `results/notebook_aligned_recovery/condensate_stage/`;
- `results/notebook_aligned_recovery/phase_stage/`;
- `results/notebook_aligned_recovery/pci_stage/`;
- `results/notebook_aligned_recovery/dgi_stage/`;
- `results/notebook_aligned_recovery/faraday_stage/`;
- `results/notebook_aligned_recovery/camera_stage/`;
- `results/notebook_aligned_recovery/multishot_stage/`.

They use `configs/notebook_v1_defaults.json`.

## Audit outputs

| Output | Generator |
| --- | --- |
| `results/performance_validation_v1/` | `scripts/run_performance_validation.py` |
| `results/thesis_numerical_consistency_v1/` | `scripts/audit_thesis_numerical_consistency.py` |
| `results/linear_approximation_audit/` | `scripts/audit_linear_approximation_validity.py` |
| `results/accumulated_snr_full_physics_audit/` | `scripts/audit_accumulated_snr_physics_and_code.py` |

## Provenance rule

No dissertation figure is complete without its config, numerical data and
metadata. Captions must distinguish provisional optical/detector inputs from
measured parameters, and every `kappa_F=1` Faraday result must remain labelled
uncalibrated.
