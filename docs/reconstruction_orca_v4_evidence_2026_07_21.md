# ORCA-Fusion reconstruction evidence, 21 July 2026

- **Status:** sealed historical synthetic evidence
- **Active consumers:** inverse-method and identifiability discussion, artifact
  verification and provenance review
- **Frozen contract:** `NA=0.080`, `kappa_F=1`, `sigma_r=1.4 e- rms`
- **Update trigger:** a demonstrated provenance, hash or transcription error;
  a new physical contract requires a new result family and evidence record
- **Retirement rule:** remove only when no current claim or regression depends
  on the sealed artifacts and the figure index is updated in the same change

## Scope

This report records the canonical synthetic reconstruction benchmark,
physical-observable recovery and credibility study for the Hamamatsu
ORCA-Fusion C14440-20UP detector contract. It covers shape-flexible dual-port
and dark-field Faraday inversion. It is not experimental validation, and the
absolute Faraday density scale remains tied to the illustrative uncalibrated
setting `kappa_F=1`.

## Frozen detector and numerical contract

| Quantity | Value |
| --- | ---: |
| Camera | ORCA-Fusion C14440-20UP |
| Effective magnification | `10` |
| Physical pixel | `6.5 um` |
| Object-plane pixel | `0.650 um` |
| Camera grid | `153 x 153` |
| Sampling | centred physical-pixel area integration |
| QE | `0.65` |
| Read noise | `1.4 e- rms` per pixel and readout |
| Reference fluence | `90 mW us` |
| Count scale at reference fluence | `220.5808727753 e- / I0 pixel` |
| Numerical aperture | `0.080` |
| Canonical propagation grid | `1024 x 1024` over `100 um` |
| Inverse grid | `306 x 306` over the same `100 um` |
| Faraday coefficient | `kappa_F=1`, illustrative and uncalibrated |

The camera values and analytic Jacobians pass through the same physical-pixel
sampler. The inverse grid retains the exact camera coordinates; the maximum
coordinate mismatch recorded by the study is
`6.78e-21 m`.

## Source freeze and provenance

The clean source freeze was created before the sealed long inverse runs:

- pre-publication local source checkpoint:
  `5e59406200f715f0891c007172f12bba7eeeec5c`;
- branch: `codex/reconstruction-v4-canonical`;
- the sealed morphology, curvature, credibility and observable numerical runs
  recorded clean generation states.

The provenance-relevant checkpoints are:

| Commit | Result | History status |
| --- | --- | --- |
| `f03802b` | morphology benchmark and sealed figures | pre-publication local checkpoint ID |
| `a05f5ec` | dual-port curvature range check | pre-publication local checkpoint ID |
| `62444d4` | dark-field curvature range check | pre-publication local checkpoint ID |
| `ea6fdc1` | active entry points, figure-label code and evidence documents | pre-publication local checkpoint ID |
| `235cfa8` | refreshed reconstruction credibility artifacts | pre-publication local checkpoint ID |
| `bd778c1` | then-active ORCA Figure 3.2 and Figure 5.2 regression artifacts | pre-publication local checkpoint ID |
| `cda577a` | production observable extraction and replay workflow | pre-publication local checkpoint ID |
| `e7e8b1d` | sealed held-out observable maps, tables and figures | pre-publication local checkpoint ID |

The prepared publication history does not retain any pre-rewrite feature
commit. The identifiers in the table, together with `5e59406`, `64a07f7` and
`df00b42`, belong to the pre-publication local DAG retained in the offline safety
refs and bundle. They remain historical provenance values but are not expected
to resolve on the eventual public branch. Public verification therefore uses
the retained config, source and artifact SHA-256 values and the deterministic
run-ID chain. Commit `498166f5d712c0a1f035fb7dd1cacd62619dc62d` is the
retrospective publication source snapshot, not the original generation commit.

The retained source-file hashes predate the repository LF contract and were
captured from a Windows CRLF worktree. Validation still requires an exact raw
hash first; its fallback accepts only the same bytes under deterministic CRLF
line-ending reconstruction. Any non-EOL source change continues to fail the
sealed-source check.

The credibility run was generated from pre-publication local checkpoint
`ea6fdc1b98c37b9f79e1bd4ef25ef058608432b1`. Its metadata records Python
3.12.13, NumPy 2.3.5, SciPy 1.18.0 and Matplotlib 3.11.1. All 13 artifact
hashes listed in that metadata were rechecked after generation and matched.

## Reduced-grid convergence

The convergence gate compares the canonical and reduced propagations on all ten
calibration and held-out analytic morphology maps. The largest recorded
relative signal-map error is `0.00374`; the largest peak error is `0.00429`.
Both are below the declared 1% limit. The earlier `204 x 204` approximation is
not used for the ORCA contract because its dark-field discrepancy exceeded the
gate.

## Frozen inverse selection

The calibration split contains five morphologies and one fixed noise
realisation per morphology. The held-out split contains five different
morphologies and two fixed realisations per morphology. Both readouts select:

```text
resolution_matched_17x5__curvature_30_um2
```

This is an 85-coefficient non-negative bilinear density with a physically
normalised thin-plate curvature weight of `30 um2`.

Stronger-weight range checks reuse the same calibration observations:

| Weight (`um2`) | Dual-port median supported-band error | Dark-field median supported-band error |
| ---: | ---: | ---: |
| 30 | 0.11439 | 0.18524 |
| 100 | 0.11725 | 0.22006 |
| 300 | 0.15744 | 0.28389 |
| 1000 | 0.25486 | 0.35880 |

Every stronger tested weight increases the median calibration error for both
readouts. The selected value is therefore bracketed against stronger
over-smoothing. This is a calibration-set check, not an additional held-out
selection.

## Held-out morphology recovery

All 60 held-out fits completed successfully: five morphologies, two draws,
three fluences and two readouts. The frozen inverse was not reselected on these
data.

| Readout | Fluence (`mW us`) | Median supported-band error | Median absolute integrated-density error |
| --- | ---: | ---: | ---: |
| Dual-port | 30 | 0.13627 | 0.05645 |
| Dual-port | 90 | 0.11178 | 0.04144 |
| Dual-port | 150 | 0.10059 | 0.04415 |
| Dark-field | 30 | 0.32608 | 0.31930 |
| Dark-field | 90 | 0.20703 | 0.18352 |
| Dark-field | 150 | 0.16346 | 0.16154 |

For the displayed held-out three-peak example at `F=90 mW us`, the
supported-band errors are `0.09097` for dual-port and `0.30528` for dark-field.
The dual-port image retains the three supported maxima more faithfully. The
dark-field reconstruction is broader and has larger background and amplitude
errors.

These errors use synthetic truth and therefore evaluate the development
contract. They are not quantities available for an unknown experimental cloud.

## Held-out physical-observable recovery

The observable study replays only the frozen 60 held-out fits; it does not
repeat calibration selection or alter the sealed morphology benchmark. Every
seed, data-Jacobian rank, parameter count and eight stored floating-point
metrics agree with the source rows. The resulting `600/600` legacy metric
checks pass.

Moments use physical cell-area weights on the common object-space density-model
support `|y| <= 27 um`, `|z| <= 7.5 um`. This support contains 7425 cells and
has area `7.9296424452e-10 m2`. It is not the camera-space fitting ROI. All 60
truth and reconstruction maps support the centroid, covariance and width
calculations. All 60 truth/reconstruction pairs also exceed the declared
fractional-anisotropy angle gate of `0.05`.

The table gives medians over ten held-out trials per readout and fluence. The
integrated-response and covariance entries are absolute relative errors; the
centroid entry is the two-dimensional position error.

| Readout | Fluence (`mW us`) | Integrated response | Centroid (`um`) | Covariance tensor |
| --- | ---: | ---: | ---: | ---: |
| Dual-port | 30 | 0.05645 | 0.1953 | 0.13866 |
| Dual-port | 90 | 0.04144 | 0.1522 | 0.06466 |
| Dual-port | 150 | 0.04415 | 0.1164 | 0.05196 |
| Dark-field | 30 | 0.31930 | 1.4331 | 0.44673 |
| Dark-field | 90 | 0.18352 | 0.4735 | 0.29081 |
| Dark-field | 150 | 0.16154 | 0.2791 | 0.31928 |

The ordered-width and aspect-ratio entries are absolute relative errors. The
angle error is the absolute axial difference modulo `pi`, evaluated only for
supported pairs.

| Readout | Fluence (`mW us`) | Major rms width | Minor rms width | Aspect ratio | Axis angle (`deg`) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Dual-port | 30 | 0.06612 | 0.15457 | 0.09187 | 0.378 |
| Dual-port | 90 | 0.03089 | 0.12712 | 0.09410 | 0.355 |
| Dual-port | 150 | 0.02452 | 0.12427 | 0.09045 | 0.242 |
| Dark-field | 30 | 0.20141 | 0.33958 | 0.12644 | 1.062 |
| Dark-field | 90 | 0.13477 | 0.22096 | 0.07491 | 0.895 |
| Dark-field | 150 | 0.14766 | 0.19144 | 0.04839 | 0.967 |

Dual-port has the smaller median integrated-response, centroid, covariance,
major-width, minor-width and angle errors at every tested fluence. Dark-field
has a smaller aspect-ratio error at `90` and `150 mW us`; this occurs while its
two width and covariance errors remain larger, so correlated size biases partly
cancel in the ratio. It is not evidence that dark-field recovers the complete
cloud size or shape more accurately.

The plotted bars are the observed minimum-to-maximum range of this finite
development ensemble, not confidence intervals. These truth-based results
evaluate the synthetic contract only. A positive reconstructed map does not by
itself establish cloud presence, and `kappa_F=1` does not make integrated
response an absolute condensate population. Relative sequence use additionally
requires stable response, detector gain and optical transfer.

After provenance paths were normalised to repository-relative form and the
retained artifact bytes were resealed under the byte-preserving result policy, the
observable numerical run has deterministic ID
`81f7c045d0745d266dbcf73cb1e908a1c016aff12b2dc53df155b6815540635b`.
This provenance-only reseal did not change the numerical arrays or tabulated
values. The run derives from morphology run
`1fb762ca7211b8df64123a346f92a21457f988ae81a5e6349a65a3249ff755e2`.
It was generated from clean pre-publication local checkpoint
`cda577ad934a6d81c03865c5c7aa660e01282401` in 1715.90 seconds using Python
3.12.13, NumPy 2.3.5, SciPy 1.18.0 and Matplotlib 3.11.1.

## Truth-independent credibility controls

The credibility study repeats the frozen inverse at the held-out three-peak
reference condition. It uses 16 conditional bootstrap draws, five curvature
weights and 12 atom-free draws per readout.

### Data support versus regularisation support

| Quantity | Dual-port | Dark-field |
| --- | ---: | ---: |
| Active coefficients | 45 | 68 |
| Effective data degrees of freedom | 34.17 | 22.88 |
| Effective prior degrees of freedom | 10.83 | 45.12 |
| Median data-mode fraction | 0.837 | 0.188 |
| Weakest data-mode fraction | 0.184 | `6.82e-5` |
| Modes with data fraction above 0.5 | 37/45 | 22/68 |
| Modes with data fraction below 0.1 | 0/45 | 27/68 |

The combined regularised systems are locally full rank. That does not mean that
all modes are measured. Most dual-port directions are data-dominated at this
condition; most dark-field directions receive their dominant support from the
curvature prior.

### Raw-channel residuals

The standardised residual rms values are `0.993` and `0.983` for the dual-port
`H` and `V` channels and `0.987` for dark-field. Residual means are within
`0.020` of zero, nearest-neighbour correlation magnitudes are below `0.025`,
and the dual-port cross-correlation is `0.0227`.

The residuals are consistent with the declared noise model. Because the same
synthetic optical operator generated and fitted the observations, this is an
internal closure check and cannot exclude unmodelled experimental structure.

### Conditional detector-noise spread

All 16 bootstrap refits completed for each readout. The recorded 68% intervals
are conditional on the fixed forward operator, basis, support, curvature prior
and solution basin. They exclude calibration and instrument-model uncertainty.

| Readout | Feature | Fitted value | Bootstrap mean | Conditional 68% interval |
| --- | --- | ---: | ---: | ---: |
| Dual-port | Integrated density | 50408 | 52372 | 51408--53523 |
| Dual-port | Peak density (`m-2`) | `5.686e14` | `5.323e14` | `5.026e14`--`5.531e14` |
| Dual-port | rms y (`um`) | 10.526 | 10.885 | 10.721--11.075 |
| Dual-port | rms z (`um`) | 1.972 | 2.176 | 2.142--2.219 |
| Dark-field | Integrated density | 64424 | 74064 | 71712--77368 |
| Dark-field | Peak density (`m-2`) | `4.891e14` | `4.519e14` | `4.267e14`--`4.715e14` |
| Dark-field | rms y (`um`) | 11.662 | 12.405 | 11.803--13.162 |
| Dark-field | rms z (`um`) | 2.578 | 2.800 | 2.732--2.878 |

Several fitted estimates lie outside their own conditional bootstrap interval.
This exposes estimator bias from the nonlinear response, positivity and
regularisation. The interval is a repeatability diagnostic, not a complete
confidence interval on the condensate.

### Prior sensitivity

Across curvature weights from `0` to `100 um2`, dual-port integrated density
changes by approximately -0.5% to +1.2%, while its peak changes by -12.2% to
+7.9%. Dark-field integrated density changes by -3.2% to +9.6%, and its peak by
-12.0% to +10.4%. Peak density, short-axis width and dark-field position or
morphology are consequently more prior-sensitive than the integrated
dual-port signal and long-axis size.

### Blank false positives

Every atom-free fit reports optimiser success, so convergence cannot be used as
evidence that atoms were present.

| Quantity | Dual-port | Dark-field |
| --- | ---: | ---: |
| Median false integrated density / reference | 0.0951 | 0.5719 |
| Upper 95% value in 12 draws | 0.1227 | 0.7066 |
| Median false peak density (`m-2`) | `4.74e13` | `1.35e14` |
| Maximum false peak density (`m-2`) | `6.58e13` | `2.08e14` |

The representative raw-count quasi-deviance improvements over a zero-density
prediction are `3381.7` for dual-port and `1262.5` for dark-field. The largest
improvements among the 12 blanks are `38.2` and `26.2`, respectively. The
representative synthetic signals rank above this small development blank set,
but these ranks are not experimental thresholds or p-values.

The dark-field blank result is the clearest warning against over-reconstruction:
a smooth positive map can be produced from noise because positivity and the
curvature prior stabilise a weak quadratic inverse. A reported reconstruction
must therefore be accompanied by raw-channel evidence, blank comparison,
data/prior support and prior sensitivity.

## Artifacts

Canonical directories:

- `results/reconstruction_morphology_benchmark_v4_orca_fusion_m10/`
- `results/reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port/`
- `results/reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field/`
- `results/reconstruction_credibility_v2_orca_fusion_m10/`
- `results/reconstruction_observables_v1_orca_fusion_m10/`

Primary figures:

- `reconstruction_quality_vs_fluence.png`
- `representative_three_peak_reconstruction_F90.png`
- `representative_credibility_F90.png`
- `data_prior_mode_support_F90.png`
- `observable_recovery_vs_fluence.png`

Reproduction commands:

Use the explicit repository interpreter defined in `reproducibility.md`; the
commands below assume its `$projectPython` variable.

```powershell
& $projectPython scripts\generate_reconstruction_morphology_benchmark.py
& $projectPython scripts\plot_reconstruction_morphology_benchmark.py
& $projectPython scripts\run_reconstruction_curvature_range_check.py --config configs\reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port.json
& $projectPython scripts\run_reconstruction_curvature_range_check.py --config configs\reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field.json
& $projectPython scripts\run_reconstruction_credibility_study.py
& $projectPython scripts\generate_reconstruction_observable_benchmark.py --validate-only
& $projectPython scripts\generate_reconstruction_observable_benchmark.py
& $projectPython scripts\plot_reconstruction_observable_benchmark.py
```

The benchmark takes approximately 45 minutes on the recorded runtime, and the
credibility study approximately 20.5 minutes. The two curvature range checks take
approximately 5 and 6.5 minutes. The held-out observable replay takes
approximately 28.6 minutes; its plot command only verifies and renders the
sealed tables. Existing artifacts should be inspected and hash-verified before
any expensive regeneration.

## Verification record

- focused reconstruction and helper tests: `128 passed`;
- complete repository suite after restoring the active ORCA result contracts:
  `259 passed in 65.97 s`;
- credibility artifact hash check: `13/13` matched;
- the then-active Figure 3.2 and Figure 5.2 numeric and PNG outputs matched the
  preceding ORCA generation byte-for-byte; only clean-branch provenance and
  vector/PDF generation metadata differed. Figure 3.2 was later retired from
  the maintained dissertation set.
- observable extraction and artifact tests: `18 passed` after adding the
  dissertation-facing numerical regression;
- complete repository suite after the observable implementation and before the
  final evidence documentation: `276 passed`; the final documented suite is
  recorded in `CODEX_HANDOFF_CURRENT.md`;
- observable replay: `60/60` fits, `600/600` legacy metric checks, six complete
  readout-by-fluence aggregate rows and zero hash failures.
