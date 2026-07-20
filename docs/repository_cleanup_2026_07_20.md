# Repository cleanup report — 20 July 2026

## Scope

Cleanup was performed on branch `work/full-performance-validation`. No merge,
checkout or write to `main` was made. The purpose was to replace the
development-time collection of trial scripts and reports with one active
dissertation and commissioning workflow.

## Active numerical contract

The repository now has one dissertation screening contract:

- `|Delta|/2pi = 1.5 GHz` reference;
- `P=1.0 mW`, `tau=90 us`, `F=90 mW us`;
- provisional `NA=0.080` and `M=4`;
- DCC3260M physical pixel `5.86 um`, object-plane pixel `1.465 um`;
- provisional `QE=0.60` and read noise `3 e-` rms;
- `2 E_rec` per scattering cycle;
- 30% condensate-depletion threshold with strict integer stopping;
- `kappa_F=1` retained only as an uncalibrated Faraday placeholder.

The contract is recorded in:

- `configs/dissertation_v2_dcc3260m.json`;
- `configs/dissertation_plots_v1.json`;
- `docs/simulation_reference_parameters.md`;
- `configs/performance_validation_v1.json`.

`configs/notebook_v1_defaults.json` remains a historical regression contract
and is no longer presented as the active detector configuration.

## Archived development branches

Fifty files or directories were moved out of the active repository to:

```text
C:\Users\jackz\Documents\non-distructive-image-method\archive\repo_cleanup_2026-07-20
```

The archive contains:

1. **Reconstruction exploration** — dual-port and dark-field Thomas-Fermi
   inversion code, configs, tests, reports and ensembles. This was a complete
   simulation experiment but was removed from the dissertation mainline.
2. **Abandoned validation and stress tests** — the deleted Section 4.6 figure
   candidates and condensate-morphology atlas.
3. **Trial and preview outputs** — M=10/M=4 trial directories, browser previews
   and cropped inspection images.
4. **Superseded transition reports** — the earlier 60 mW us, camera-selection,
   recoil-update and July 14 placeholder-verification records.
5. **Superseded planning documents** — obsolete branch manifests,
   dissertation structures, optimisation-readiness and calibration-roadmap
   notes.
6. **Superseded numerical reports** — long hand-written reports that mixed the
   old detector/sampling contract with current results.

The archive is intentionally outside the Git repository. Git history also
retains all previously tracked documents.

## Active repository surface

The current active documentation is limited to:

- model architecture and physics;
- the current screening contract;
- figure and result provenance;
- canonical and numerical-consistency audits;
- figure-language conventions;
- the detailed laboratory measurement plan;
- reproducibility and release instructions.

The active dissertation result tree contains:

- Figure 3.2;
- Figure 4.2;
- Figure 5.1;
- the Figure 5.2 scan, dual-port heatmap and operating-band view;
- Figure 5.4;
- the full four-mode multishot comparison.

Trial figures and abandoned reconstruction outputs are no longer mixed into
`results/dissertation_plots_v1/`.

## Consistency problems found and fixed

### Two detector contracts

The canonical performance config still used `QE=0.5851647` and
`sigma_r=7 e-`, while the dissertation configs used `QE=0.60` and
`sigma_r=3 e-`. The gate and tests now use the active provisional
`0.60/3 e-` contract.

### Stale M=2 multishot data

The active config specified `M=4`, but the stored full-multishot output still
used a 34x34 camera grid, bin size 30 and a 105-pixel ROI. All retained outputs
were regenerated from the active config. The current result uses:

- 68x68 camera sampling;
- bin size 15;
- a common 228-pixel matched ROI.

### Broken clean regeneration

After Figure 3.2 was regenerated, the full-multishot generator still requested
the obsolete CSV column `SNR_total`. The active column is `SNR_acc`. The
consumer was corrected, so the retained result chain now regenerates from
scratch rather than only working with an old CSV.

### Incomplete reproduction entry point

`run_all_dissertation_figures.py` now includes Figures 4.2, 5.1, 5.2 and 5.4
and finishes with the canonical gate. The old placeholder Faraday optimisation
plot is retained as regression material but is no longer regenerated as an
approved dissertation output.

## Current canonical result

At the reference operating point, 10 frames are accepted. The loss after the
tenth pulse is `0.2806023`; the next pulse would give `0.3082904`.

Accumulated matched-template SNR over the evolving sequence:

| Mode | Shot only | Shot + read |
| --- | ---: | ---: |
| PCI | 104.0982 | 103.6687 |
| DGI | 54.3620 | 34.3033 |
| Faraday dark-field | 55.7783 | 34.7088 |
| Faraday dual-port | 107.4820 | 106.5555 |

These are 228-pixel matched-template results. They must not be compared as if
they were the single-frame central-`3x3` SNR used in the Chapter 5 camera
figures. Faraday values remain uncalibrated.

## Validation completed

- canonical gate: passed;
- focused current-contract regression set: 16 passed;
- full suite after cleanup and current-output regeneration: 139 passed;
- reproduction dry run: passed and lists every current figure plus the
  canonical gate.

## Deliberately retained legacy material

The original notebook, exported notebook sections, regression baselines,
notebook-aligned recovery outputs and legacy numerical-correction audit remain
in the repository. They are retained because they document provenance and
protect the migrated equations. They are explicitly separated from the active
screening contract.

This cleanup is published as part of the Dissertation v2 update.
