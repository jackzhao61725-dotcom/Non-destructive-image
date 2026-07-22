# Figure and data index

- **Status:** single result-to-generator-to-config status index
- **Active consumers:** dissertation figure review, reproduction and repository
  retention decisions
- **Update trigger:** a result family, generator, config, dissertation consumer
  or scientific-status label changes
- **Retirement rule:** replace only with one indexed successor that covers every
  retained result family

This is the single index from a maintained output to its generator, config,
quantity and scientific status. Each generator writes numerical data and
metadata alongside the rendered figure. A stored output is not current merely
because it is reproducible.

## Dissertation figures

The Figure 4.2-5.4 family uses signed `kappa_F=-45/91`, effective
`NA=0.130`, the ORCA-Fusion Ultra quiet manufacturer-typical
`sigma_r=0.7 e- rms` scenario, a `90-300 mW us` fluence scan, representative
conditions `F=90, 150, 300 mW us`, and `F=300 mW us` as the first-frame
reference. These are screening inputs, not installed-apparatus measurements.

| Output | Generator | Config | Quantity and scope | Status |
| --- | --- | --- | --- | --- |
| `results/dissertation_plots_v2_orca_fusion/figure_4_2/figure_4_2.{png,pdf,svg}` | `scripts/generate_figure_4_2.py` | `configs/figure_4_2.json` | noiseless finite-aperture PCI, DGI and signed Faraday outputs at effective `NA=0.130` | active regenerated figure |
| `results/dissertation_plots_v2_orca_fusion/figure_5_1/figure_5_1.{png,pdf,svg}` | `scripts/generate_figure_5_1.py` | `configs/figure_5_1.json` | analytic `SNR_5x5` over `F=90-300 mW us` and noisy frames at `90, 150, 300 mW us` | active screening figure; `F=300 mW us` is the first-frame reference |
| `results/dissertation_plots_v2_orca_fusion/figure_5_2/figure_5_2_dual_port_heatmap.{png,pdf,svg}` | `scripts/generate_figure_5_2.py` | `configs/figure_5_2.json` | strict usable-frame screen over detuning and `F=90-300 mW us` | frozen Version 1 sequence screen; canonical regeneration blocked pending heating replacement |
| `results/dissertation_plots_v2_orca_fusion/figure_5_2/figure_5_2_operating_band.{png,pdf,svg}` | `scripts/generate_figure_5_2.py` | `configs/figure_5_2.json` | fixed-1.5-GHz comparison of Version 1 `N_dep` and `N_use` | frozen Version 1 operating band; canonical regeneration blocked pending heating replacement |
| `results/dissertation_plots_v2_orca_fusion/figure_5_4/figure_5_4.{png,pdf,svg}` | `scripts/generate_figure_5_4_snr_panel.py` | `configs/figure_5_4.json` | 15-state dual-port sequences at `F=90, 150, 300 mW us` | frozen Version 1 illustration; canonical regeneration blocked; not observable-recoverability evidence |

The `figure_5_2` directory also contains a combined overview figure and the
complete scan CSV/JSON. Which rendered panel receives the final dissertation
figure number is a writing decision; the numerical source remains one scan.

## Reconstruction figures

Every reconstruction output in this section is sealed historical method
evidence under `NA=0.080`, `kappa_F=1` and `sigma_r=1.4 e- rms`. It is retained
without regeneration and must not be cited as quantitative performance under
the active screening contract.

| Output | Generator | Config | Quantity and scope | Status |
| --- | --- | --- | --- | --- |
| `results/reconstruction_morphology_benchmark_v4_orca_fusion_m10/reconstruction_quality_vs_fluence.{png,pdf,svg}` | `scripts/plot_reconstruction_morphology_benchmark.py` | `configs/reconstruction_morphology_benchmark_v4_orca_fusion_m10.json` | held-out pupil-supported full-field diagnostic | historical frozen diagnostic, not a primary dissertation result |
| `results/reconstruction_morphology_benchmark_v4_orca_fusion_m10/representative_three_peak_reconstruction_F90.{png,pdf,svg}` | same | same | one truth, raw observation and latent-fit comparison | frozen diagnostic; no image-recovery claim |
| `results/reconstruction_credibility_v2_orca_fusion_m10/representative_credibility_F90.{png,pdf}` | `scripts/run_reconstruction_credibility_study.py` | `configs/reconstruction_credibility_v2_orca_fusion_m10.json` | raw observation, latent fit, conditional detector-noise spread and residuals | frozen diagnostic under the recorded operator |
| `results/reconstruction_credibility_v2_orca_fusion_m10/data_prior_mode_support_F90.{png,pdf}` | same | same | local division of inverse support between camera data and curvature prior | frozen diagnostic under the recorded operator |
| `results/reconstruction_observables_v1_orca_fusion_m10/observable_recovery_vs_fluence.{png,pdf,svg}` | `scripts/plot_reconstruction_observable_benchmark.py` | `configs/reconstruction_observables_v1_orca_fusion_m10.json` | truth-known errors for `A`, centroid and width; any extra plotted quantities are diagnostics; bars are finite-ensemble ranges | historical sealed evidence at `NA=0.080`, `kappa_F=1`, `sigma_r=1.4 e- rms` |

The representative reconstruction figures are synthetic. The first benchmark
figure may use truth error because truth is available only in the development
study. The credibility figures do not use truth as their evidence measure and
must be interpreted together with the blank and prior-sensitivity tables. The
observable figure is also a truth-based synthetic development result. Its
complete covariance errors and supported denominators remain in
`observable_recovery_summary.csv` even though the compact figure does not plot
the covariance panel.

## Supporting physics plot

| Output | Generator | Config | Scope |
| --- | --- | --- | --- |
| `results/dissertation_plots_v2_orca_fusion/detuning_tradeoff/detuning_tradeoff.svg` | `scripts/generate_detuning_tradeoff_plot.py` | `configs/dissertation_plots_v2_orca_fusion.json` | normalised dispersive-signal and scattering trends |

## Historical notebook-aligned figures

The complete `results/notebook_aligned_recovery/` family is retained as
historical regression and interpretation evidence, not as a current detector
prediction. Its generated subdirectories are:

- `condensate_stage/` and `condensate_three_view/`;
- `phase_stage/`, `pci_stage/` and `dgi_stage/`;
- `faraday_stage/` and `faraday_camera_panel/`;
- `camera_stage/` and `noisy_camera_stage/`;
- `multishot_stage/` and `noisy_multishot_filmstrip/`.

They use `configs/notebook_v1_defaults.json`.

## Result-directory roles

- `results/dissertation_plots_v2_orca_fusion/` holds detector-dependent figure
  data under the signed-`kappa_F`, effective-`NA=0.130`, Ultra quiet screening
  contract. Figure 4.2 and Figure 5.1 are current; Figures 5.2 and 5.4 are frozen
  Version 1 sequence outputs whose canonical regeneration is blocked pending
  the approved heating replacement.
- `results/reconstruction_morphology_benchmark_v4_orca_fusion_m10/` and the
  associated curvature-range directories retain the frozen inverse selection.
- `results/reconstruction_credibility_v2_orca_fusion_m10/` retains residual,
  data/prior, detector-noise, prior-sensitivity and blank controls.
- `results/reconstruction_observables_v1_orca_fusion_m10/` retains the sealed
  60-fit truth-known observable replay.
- `results/notebook_aligned_recovery/` is historical regression evidence only.

## Output admission policy

A new maintained result family requires:

1. an explicit config and deterministic generator;
2. machine-readable values and metadata with code and parameter provenance;
3. a current dissertation, regression or calibration consumer;
4. a regression check when the output supports a dissertation claim;
5. a distinct directory only when it supersedes a different frozen contract.

Do not retain trial crops, browser previews, alternative plotting experiments
or regenerable renders without a current consumer. Remove superseded outputs
from the active tree; Git history is the archive.

## Provenance rule

No dissertation figure is complete without its config, numerical data,
metadata and scientific status. Captions distinguish provisional inputs,
theoretical estimates, calibrated parameters and held-out validation. Active
screening outputs use signed `kappa_F=-45/91`, effective `NA=0.130` and the
manufacturer-typical Ultra quiet `0.7 e- rms` scenario. Every stored
`NA=0.080`, `kappa_F=1`, `sigma_r=1.4 e- rms` reconstruction result is labelled
as historical frozen evidence, not merely as uncalibrated.
