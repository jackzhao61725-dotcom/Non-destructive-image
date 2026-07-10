# Figure Index

This index lists approved dissertation-facing figures and plots available on
`main` at the time the reproducibility layer was prepared.

## Notebook-Aligned Recovery Figures

| Figure | Generating script | Config | Type | Calibration status | Metadata | Use |
| --- | --- | --- | --- | --- | --- | --- |
| `results/notebook_aligned_recovery/condensate_stage/condensate_density_stage.svg` | `scripts/recover_notebook_condensate_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/condensate_stage/metadata.json` | main text |
| `results/notebook_aligned_recovery/phase_stage/scalar_phase_stage.svg` | `scripts/recover_notebook_phase_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/phase_stage/metadata.json` | main text |
| `results/notebook_aligned_recovery/pci_stage/pci_image_stage.svg` | `scripts/recover_notebook_pci_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/pci_stage/metadata.json` | main text or appendix |
| `results/notebook_aligned_recovery/dgi_stage/dgi_image_stage.svg` | `scripts/recover_notebook_dgi_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/dgi_stage/metadata.json` | main text or appendix |
| `results/notebook_aligned_recovery/faraday_stage/faraday_dark_field_stage.svg` | `scripts/recover_notebook_faraday_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1, `kappa_F = 1.0` | `results/notebook_aligned_recovery/faraday_stage/metadata.json` | main text |
| `results/notebook_aligned_recovery/faraday_stage/faraday_dual_port_signal_stage.svg` | `scripts/recover_notebook_faraday_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1, `kappa_F = 1.0` | `results/notebook_aligned_recovery/faraday_stage/metadata.json` | main text |
| `results/notebook_aligned_recovery/camera_stage/camera_deterministic_stage.svg` | `scripts/recover_notebook_camera_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/camera_stage/metadata.json` | methods or appendix |
| `results/notebook_aligned_recovery/noisy_camera_stage/noisy_camera_stage.svg` | `scripts/recover_notebook_noisy_camera_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1, explicit RNG | `results/notebook_aligned_recovery/noisy_camera_stage/metadata.json` | methods or appendix |
| `results/notebook_aligned_recovery/multishot_stage/multishot_sequence_stage.svg` | `scripts/recover_notebook_multishot_stage.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1 | `results/notebook_aligned_recovery/multishot_stage/metadata.json` | main text |
| `results/notebook_aligned_recovery/noisy_multishot_filmstrip/noisy_multishot_pci_filmstrip.svg` | `scripts/recover_notebook_noisy_multishot_filmstrip.py` | `configs/notebook_v1_defaults.json` | notebook-aligned recovery | uncalibrated V1, explicit RNG | `results/notebook_aligned_recovery/noisy_multishot_filmstrip/metadata.json` | appendix unless visually approved |

## Model Extension Figures

| Figure | Generating script | Config | Type | Calibration status | Metadata | Use |
| --- | --- | --- | --- | --- | --- | --- |
| `results/notebook_aligned_recovery/condensate_three_view/condensate_three_view.svg` | `scripts/generate_condensate_three_view.py` | `configs/notebook_v1_defaults.json` | notebook-aligned condensate-model extension | uncalibrated V1 | `results/notebook_aligned_recovery/condensate_three_view/metadata.json` | appendix |

## Representative V1 Plots

| Figure | Generating script | Config | Type | Calibration status | Metadata | Use |
| --- | --- | --- | --- | --- | --- | --- |
| `results/faraday_optimisation_v1/detuning_tradeoff.svg` | `scripts/generate_dissertation_results.py` | `configs/dissertation_results_v1.json` | representative V1 plot | uncalibrated, `kappa_F = 1.0` | `results/faraday_optimisation_v1/metadata.json` | appendix or workflow demonstration |
| `results/faraday_optimisation_v1/intensity_tradeoff.svg` | `scripts/generate_dissertation_results.py` | `configs/dissertation_results_v1.json` | representative V1 plot | uncalibrated, `kappa_F = 1.0` | `results/faraday_optimisation_v1/metadata.json` | appendix or workflow demonstration |
| `results/faraday_optimisation_v1/exposure_time_tradeoff.svg` | `scripts/generate_dissertation_results.py` | `configs/dissertation_results_v1.json` | representative V1 plot | uncalibrated, `kappa_F = 1.0` | `results/faraday_optimisation_v1/metadata.json` | appendix or workflow demonstration |

## Pending Figure Workflows

| Figure | Branch | Status |
| --- | --- | --- |
| `results/dissertation_plots_v1/detuning_tradeoff/detuning_tradeoff.svg` | `work/detuning-tradeoff-plot` | pending approval and merge into `main` |

## Provenance Rule

No figure without parameter provenance. Every approved figure or plot should
have metadata documenting config files, parameters, units, normalisation, and
calibration status.
