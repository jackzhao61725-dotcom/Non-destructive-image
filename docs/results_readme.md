# Results and figure generation

Generated outputs are retained only when they support the historical regression
baseline, the current dissertation figures or the canonical numerical gate.
Every active result directory contains machine-readable metadata.

## Dissertation figures

| Result directory | Content | Generator | Config |
| --- | --- | --- | --- |
| `results/dissertation_plots_v1/accumulated_snr_invariance/` | Figure 3.2 and parameter register | `generate_accumulated_snr_invariance_plot.py` | `dissertation_plots_v1.json` |
| `results/dissertation_plots_v1/figure_4_2/` | noiseless four-readout comparison | `generate_figure_4_2.py` | `figure_4_2.json` |
| `results/dissertation_plots_v1/figure_5_1/` | single-frame Faraday SNR and camera examples | `generate_figure_5_1.py` | `figure_5_1.json` |
| `results/dissertation_plots_v1/figure_5_2/` | usable-frame screen, dual-port map and operating band | `generate_figure_5_2.py` | `figure_5_2.json` |
| `results/dissertation_plots_v1/figure_5_4/` | 25-frame dual-port filmstrip | `generate_figure_5_4_snr_panel.py` | `figure_5_4.json` |
| `results/dissertation_plots_v1/full_multishot_accumulated_snr/` | four-mode evolving matched-template SNR data | `generate_full_multishot_accumulated_snr.py` | `dissertation_plots_v1.json` |

The active detector/sampling contract is `M=4`, `QE=0.60`,
`sigma_r=3 e-` rms. The reference operating point is
`|Delta|/2pi=1.5 GHz`, `P=1 mW`, `tau=90 us`.

## Evidence and audits

| Directory | Purpose |
| --- | --- |
| `results/performance_validation_v1/` | canonical gate and full provenance ledger |
| `results/thesis_numerical_consistency_v1/` | current tables plus explicit legacy corrections |
| `results/linear_approximation_audit/` | weak-phase and small-angle numerical checks |
| `results/accumulated_snr_full_physics_audit/` | machine-readable comparison of clean-loss and evolving sequence assumptions |
| `results/notebook_aligned_recovery/` | historical notebook-stage regression evidence |

The `faraday_optimisation_v1` directory is a retained workflow/regression
example, not the current dissertation operating-region result.

## Focused regeneration

From the repository root:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"

python scripts\generate_accumulated_snr_invariance_plot.py --config configs\dissertation_plots_v1.json
python scripts\generate_full_multishot_accumulated_snr.py --config configs\dissertation_plots_v1.json
python scripts\generate_figure_4_2.py --config configs\figure_4_2.json
python scripts\generate_figure_5_1.py --config configs\figure_5_1.json
python scripts\generate_figure_5_2.py --config configs\figure_5_2.json
python scripts\generate_figure_5_4_snr_panel.py --config configs\figure_5_4.json
python scripts\audit_thesis_numerical_consistency.py --config configs\thesis_numerical_contract_v1.json
python scripts\run_performance_validation.py --skip-prerequisites
```

For a full ordered run, use
`python scripts\run_all_dissertation_figures.py`.

## Output policy

Do not add trial crops, browser previews or alternative magnification trials to
the maintained result tree. A new result family needs:

1. an explicit config;
2. a deterministic generator;
3. metadata with parameter and calibration provenance;
4. a regression check if the result is dissertation-facing;
5. a distinct output directory when it supersedes rather than reproduces an
   existing contract.
