# Canonical Faraday Full-Multishot Matched-ROI SNR

## Scope

This result extends the existing canonical PCI/DGI full-multishot pipeline to
the two Version 1 Faraday readouts:

- crossed-analyser dark-field intensity;
- dual-port difference of the two linear-polarisation outputs.

It does not use `evaluate_faraday_operating_point(...)` or the representative
Faraday optimisation config. Every mode uses the same condensate, optical
point, camera scale, fixed ROI, heating/reabsorption sequence, and strict
integer stopping rule in
`scripts/generate_full_multishot_accumulated_snr.py`.

## Canonical Parameter Register

| Quantity | Value |
| --- | ---: |
| Initial condensate atom number | `25000` |
| Initial peak x-axis column density | `5.3759624525784675e14 m^-2` |
| Saturation intensity | `597.9627307559916 W m^-2` |
| `|Delta|/2pi` | `1.5 GHz` |
| Probe power | `3.5 mW` |
| Exposure time | `40 us` |
| Imaging axis | `x` |
| Quantum efficiency | `0.4` |
| Read noise | `7 e- rms` per port and camera pixel |
| Matched ROI | `228` binned camera pixels |
| Depletion threshold | `30%` of initial condensate number |
| Sequence model | heating plus initial-density reabsorption |
| Stopping semantics | strict integer; threshold-crossing pulse excluded |
| Faraday factor | `kappa_F = 1.0` placeholder |

At this point the continuous heating/reabsorption crossing is `13.7576`
pulses. The pipeline accepts 13 frames: depletion after frame 13 is
`0.2836815`; one more pulse would give `0.3052166` and is excluded.

## Matched-ROI Observables

For dark-field Faraday, the electron template is the crossed-analyser count
map relative to the ideal zero background. Per-pixel variances are

```text
shot only:  N_dark
shot+read:  N_dark + sigma_read^2
```

For dual-port Faraday, the electron template is the count difference
`N_V-N_U`. Independent Poisson ports and independent detector readouts give

```text
shot only:  N_U + N_V
shot+read:  N_U + N_V + 2 sigma_read^2
```

This is algebraically the notebook's normalised-difference error propagation:
the common `N_U+N_V` denominator cancels from the per-pixel SNR. The full result
uses the diagonal-covariance matched-template convention
`sqrt(sum_i signal_i^2 / variance_i)` over the same fixed 228-pixel ROI used
for PCI and DGI, followed by framewise RMS accumulation.

## Canonical Results

| Mode | Noise model | Accumulated matched-ROI SNR |
| --- | --- | ---: |
| Faraday dark-field | shot-noise only | `64.46872450095157` |
| Faraday dark-field | shot plus read | `22.163032104041168` |
| Faraday dual-port | shot-noise only | `124.23007404156489` |
| Faraday dual-port | shot plus read | `118.89881207407808` |

The dark-field channel is strongly read-noise limited because the crossed
analyser rejects the bright carrier. The dual-port channel retains the photon
flux in two bright ports; read noise therefore produces a much smaller
reduction. These mechanisms are consistent with the corresponding DGI and PCI
behaviour under the same sequence.

## Calibration Boundary

These four values are **Version 1 representative and uncalibrated**. Their
absolute scale follows the phenomenological convention

```text
theta_F = kappa_F * phi,  kappa_F = 1.0.
```

They complete the structural four-mode comparison but are not calibrated
Faraday predictions. They must not be presented as an experimentally validated
operating point until `kappa_F` and the relevant vector-polarisability response
have been fitted to laboratory data.

Additional omissions are analyser leakage/extinction, correlated technical
noise, reference-image noise, inter-port gain imbalance, and density-updated
reabsorption within the sequence.

## Reproduction

```powershell
$env:PYTHONPATH = "src;."
$env:PYTHONUTF8 = "1"
python scripts\generate_full_multishot_accumulated_snr.py --config configs\dissertation_plots_v1.json
```

The machine-readable parameter ledger is
`results/dissertation_plots_v1/full_multishot_accumulated_snr/faraday_canonical_reference_at_1p5GHz.csv`.
