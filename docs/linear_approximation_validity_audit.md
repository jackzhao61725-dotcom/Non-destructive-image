# Linear Approximation Validity Audit

This audit checks whether the current Version 1 dissertation simulation and interpretation are justified in using small-angle or far-detuned language. It is an audit only: no simulator physics, helper API, notebook section, baseline, or existing recovery output is changed.

## Parameter Provenance

- Canonical config: `configs\notebook_v1_defaults.json`
- Calibration status: audit only; no recalibration performed
- Detuning convention: `delta = 2 * detuning_hz * 2*pi / Gamma_rad_per_s`
- Git commit audited: `d6b0ac23df4e0a492ea407c5d1ebdd7ce94a7ed0`
- Detuning plot config present on this branch: `True`

## Approximation Sites

| Site | Topic | Code status | Approximation role |
| --- | --- | --- | --- |
| `notebook_sections/03_light_atom_interaction.py cells 9-12` | scalar phase and far-detuned scaling | exact lineshape for phi_peak and residual OD | interpretation and regime classification |
| `notebook_sections/04_pci.py cells 15-18 and 22` | PCI transfer and SNR | exact exp(i phi) propagation for images; linear formula only for ideal SNR/reference tangent | interpretation and idealised SNR comparison |
| `src/non_destructive_image/imaging.py simulate_pci_image and simulate_dgi_image` | scalar field propagation | exact complex phase field | not applicable |
| `notebook_sections/06_faraday.py cells 84-89` | Faraday rotation readout | field recombination with exp(+/- i theta_F), then exact intensities/ratio | interpretation of dark-field and dual-port response |
| `src/non_destructive_image/imaging.py simulate_faraday_image` | Faraday field propagation | exact finite-rotation field expression for Version 1 phenomenological theta_F | not applicable to propagation; theta_F model itself is phenomenological |
| `src/non_destructive_image/analysis.py evaluate_faraday_operating_point` | optimisation signal proxy | scalar proxy proportional to theta_F, not an image-field simulation | single-variable optimisation proxy, not a final camera/image prediction |
| `notebook_sections/09_multishot_simulation.py cells 40-46` | multishot signal/loss interpretation | exact phase lineshape, exact PCI transfer for realistic SNR, and exact scattering denominator | physical explanation and idealised invariance arguments |

## Canonical Phase And Rotation Ranges

- max |phi| = 0.202942 rad
- central |phi| = 0.202942 rad
- 95th percentile |phi| over nonzero cloud pixels = 0.188189 rad
- 99th percentile |phi| over nonzero cloud pixels = 0.200082 rad
- max |theta_F| = 0.202942 rad
- central |theta_F| = 0.202942 rad
- 95th percentile |theta_F| over nonzero cloud pixels = 0.188189 rad
- 99th percentile |theta_F| over nonzero cloud pixels = 0.200082 rad

Because the current Version 1 Faraday convention uses `kappa_F = 1.0`, the canonical theta_F range is identical to the scalar phase range.

## Small-Angle Error Summary

- At peak |phi|, `exp(i phi)` versus `1 + i phi` has relative field error 2.0569%.
- At peak |theta_F|, `sin(theta_F)` versus `theta_F` has relative error 0.6897%.
- At peak |theta_F|, dark-field `sin^2(theta_F)` versus `theta_F^2` has relative error 1.3842%.
- At peak |theta_F|, dual-port `sin(2 theta_F)` versus `2 theta_F` has relative error 2.7994%.

These errors are small enough for qualitative scaling language at the canonical operating point, but they are not zero. Quantitative image formation should continue to use the exact expressions already present in the recovered code.

## Far-Detuned Scaling

At 1.5 GHz, delta = 101.695. The relative error of replacing the exact phase response `delta/(1+delta^2)` with `1/delta` is 0.0097%. The relative error of replacing residual OD with `1/delta^2` is 0.0097%. The scattering scaling error, including the audited saturation parameter, is 0.0099%.

The 1/Delta and 1/Delta^2 statements are therefore acceptable as far-detuned scaling trends over the audited detuning range, not as exact equalities.

## Validity Classification

| Topic | Classification | Reason |
| --- | --- | --- |
| scalar phase propagation | Not applicable because code already uses exact expression | Images use exp(i phi). The peak first-order field-expansion error is 2.057%, so first-order propagation should not replace the exact field. |
| PCI interpretation | Safe to marginal as an interpretive tangent only | Canonical peak |phi|=0.203 rad and p99=0.200 rad are below the notebook phi<0.5 PCI-linear guide, but finite-phase curvature is measurable. |
| Faraday rotation angle | Marginal but acceptable with caveat | kappa_F=1 gives peak |theta_F|=0.203 rad. This is weak rotation, not infinitesimal rotation. |
| Faraday dark-field signal | Marginal but acceptable with caveat | At peak theta_F, sin^2(theta) vs theta^2 differs by 1.384%. Use exact expression in simulation and quadratic only as scaling. |
| Faraday dual-port signal | Marginal but acceptable with caveat | At peak theta_F, sin(2 theta) vs 2 theta differs by 2.799%; the exact dual-port expression should be cited for quantitative results. |
| detuning scaling plot | Safe far-detuned scaling at >=0.5 GHz for trend statements | For the audited detunings, |delta| is large and 1/Delta, 1/Delta^2 errors are small; wording should still call them asymptotic trends. |
| multishot signal/loss interpretation | Not applicable because code already uses exact expression | Sequence recovery uses exact scattering denominator, exact phase expression, and full PCI transfer for realistic SNR; linear invariance arguments are explanatory. |

## Dissertation Wording Recommendations

- Say: the simulation propagates the full complex field `exp(i phi)`; linear phase response should be used only as an interpretive scaling statement.
- Say: the canonical operating point has finite phase/rotation values, so exact numerical expressions are retained for images and quantitative comparisons.
- Say: Faraday rotation remains a weak-rotation regime for Version 1, but not a strictly infinitesimal-rotation limit.
- Say: dark-field Faraday scales approximately as theta_F squared and dual-port Faraday approximately as 2 theta_F only in the small-angle interpretation; the plotted/recovered fields use exact expressions.
- Say: 1/Delta and 1/Delta^2 are far-detuned trends, not exact equalities across all plotted detunings.
- Avoid claiming that representative V1 optimisation outputs are final calibrated operating-point predictions; no absorption/RAI calibration or kappa_F calibration has been applied.

## Output Files

- `results/linear_approximation_audit/linear_approximation_summary.json`
- `results/linear_approximation_audit/phase_rotation_ranges.csv`
- `results/linear_approximation_audit/small_angle_error_table.csv`
- `results/linear_approximation_audit/detuning_scaling_error_table.csv`
- `results/linear_approximation_audit/metadata.json`
