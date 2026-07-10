# Notebook-Aligned Parameter Consistency Audit

Audited source state: `main` at `fcd9865` (`Recover notebook noisy multishot filmstrip`).

This audit checks whether the notebook-aligned recovery chain uses `configs/notebook_v1_defaults.json` as the canonical parameter source, and whether recovery scripts, metadata, summaries, and tests remain consistent with that source.

The full parameter inventory is saved at:

- `results/notebook_aligned_recovery/parameter_inventory.csv`

## Scope

Inspected recovery inputs and outputs:

- `configs/notebook_v1_defaults.json`
- `scripts/recover_notebook_condensate_stage.py`
- `scripts/recover_notebook_phase_stage.py`
- `scripts/recover_notebook_pci_stage.py`
- `scripts/recover_notebook_dgi_stage.py`
- `scripts/recover_notebook_faraday_stage.py`
- `scripts/recover_notebook_camera_stage.py`
- `scripts/recover_notebook_noisy_camera_stage.py`
- `scripts/recover_notebook_multishot_stage.py`
- `scripts/recover_notebook_noisy_multishot_filmstrip.py`
- notebook-aligned recovery metadata and summary JSON files under `results/notebook_aligned_recovery/`
- notebook recovery regression tests under `tests/regression/`
- `docs/notebook_aligned_recovery_status.md`

No simulator code, notebook sections, baselines, or existing recovery result data were changed.

## Main Findings

### 1. Single source of notebook defaults

`configs/notebook_v1_defaults.json` is the effective single source for the active notebook-aligned recovery workflow.

The active recovery scripts accept a `--config` argument and load this file before computing stage outputs. The main physical, optical, camera, stochastic, and multishot parameters are centralized there:

- condensate parameters: atom number, mass, scattering length, trap frequencies, temperature
- grid parameters: grid size and field of view
- light-atom parameters: wavelength, linewidth, resonant cross section, scalar phase detuning
- imaging parameters: optical geometry, PCI phase plate, DGI optical depth, Faraday `kappa_F`
- camera parameters: power, exposure, bin size, quantum efficiency, read noise
- stochastic parameters: explicit RNG seeds
- multishot parameters: pulse duration, power, loss threshold, maximum shots, collision-loss coefficient
- filmstrip display parameters: selected frame indices and display colour limits

### 2. Scripts generally load from config

The active stage recovery scripts load notebook defaults from the config and then compute derived quantities from those values. This is the right structure for traceability: changing a notebook default in one place should propagate through the recovery chain.

Hard-coded values in the scripts fall into three categories:

- **Display and file-output constants:** figure size, line width, output folder names, CSV column names, colour labels.
- **Notebook formula constants:** algebraic coefficients such as `15`, `8*pi/15`, `4/3`, `2*pi`, and the Bose critical-temperature prefactor `0.94`.
- **Numerical procedure constants:** root-solver brackets such as `[3e4, 1e7]` atoms for the self-consistent total-atom calculation.

Most hard-coded values are not independent physical defaults. The two values most worth documenting or moving into config later are:

- `0.94`, the Bose critical-temperature prefactor used in the notebook formula.
- `[3e4, 1e7]`, the self-consistent total-atom root-solver bracket.

These are not currently causing numerical inconsistency, but they are traceability weak points.

### 3. Metadata and summaries match config

Stage metadata and summary files are consistent with the current config for the parameters they record:

- condensate summary records `N0 = 25000`, `166Er`, `a_s = 72 a0`, trap frequencies `[293, 14, 233] Hz`
- phase summary records detuning `1.5e9 Hz`
- PCI summary records `t_p = 0.95` and `theta = pi/2`
- DGI summary records `OD = 4.0`
- Faraday summary records `kappa_F = 1.0`
- camera summary records `probe_power = 2.0 mW`, `exposure = 100 us`, `QE = 0.4`, `read_noise = 7 e-`, `bin_size = 15`
- noisy camera summary records explicit seed `7`
- multishot summary records `power = 3.5 mW`, `pulse = 40 us`, `max_shots = 400`, `loss_limit = 0.3`, `eta_coll = 1.3`
- noisy filmstrip summary records selected frames `[0, 5, 10, 14]` and explicit seed `7`

The single-frame camera and multishot stages intentionally use different probe powers and exposure/pulse durations. This is a stage-specific notebook distinction, not a mismatch.

### 4. Tests are stable but duplicate reference numbers

The notebook recovery tests use expected numerical values derived from the current recovery outputs/config. This is appropriate for regression testing, but it means the tests intentionally duplicate some config-derived results.

Classification:

- **Stable regression expectation:** acceptable.
- **Not a second source of scientific defaults:** tests should be updated if `configs/notebook_v1_defaults.json` is intentionally changed.

No test was found to introduce a separate physical model or independent notebook default.

### 5. Stage-specific differences are legitimate

The following stage-specific differences are expected:

- `camera_recovery.probe_power_mw = 2.0` and `camera_recovery.default_exposure_s = 1e-4` are used for single-frame deterministic/noisy camera recovery.
- `multishot_recovery.power_mw = 3.5` and `multishot_recovery.pulse_duration_s = 4e-5` are used for continuous-imaging/multishot recovery.
- DGI uses `OD = 4.0` as a DGI-only reference attenuation parameter.
- Faraday uses `kappa_F = 1.0` as the Version 1 phenomenological placeholder.
- Filmstrip selected frames and colour limits are display-level choices.

These are documented in the inventory as `STAGE_SPECIFIC`, not inconsistent.

## Parameter Traceability Assessment

| Parameter group | Status | Notes |
| --- | --- | --- |
| Condensate defaults | CONSISTENT | Atom number, species, scattering length, trap frequencies, and temperature are centralized. |
| Grid and units | CONSISTENT | Internal coordinates are SI; plotting and CSV labels convert to micrometers where appropriate. |
| Scalar phase defaults | CONSISTENT | Detuning and atomic constants are config-driven. |
| PCI defaults | CONSISTENT | Phase-plate transmission and phase shift are config-driven. |
| DGI defaults | STAGE_SPECIFIC | `OD = 4.0` is DGI-only and documented. |
| Faraday defaults | STAGE_SPECIFIC | `kappa_F = 1.0` remains a placeholder and is documented. |
| Deterministic camera defaults | STAGE_SPECIFIC | Single-frame camera settings differ from multishot settings by notebook design. |
| Stochastic camera defaults | CONSISTENT | Explicit seeds are used for reproducible recovery; exact notebook-global RNG state remains ambiguous. |
| Multishot defaults | STAGE_SPECIFIC | Power, pulse duration, loss threshold, and frame count match notebook recovery intent. |
| Display choices | STAGE_SPECIFIC | Figure limits and colour ranges are display-only and should not be treated as physical outputs. |
| Formula/procedure constants | MISSING_FROM_CONFIG | The `0.94` Tc prefactor and solver brackets are still embedded in scripts. |

## Units Check

The recovery workflow is internally SI-based:

- length: meters internally, micrometers for plotted axes and lineout labels
- detuning: Hz in config, dimensionless detuning derived using linewidth
- power: milliwatts in config where notebook used milliwatt-scale inputs
- exposure/pulse time: seconds
- trap frequencies: Hz, converted to angular frequencies where equations require it
- read noise: electrons RMS
- density/column density: SI-derived arrays from the notebook Thomas-Fermi convention

No unit conflict was found in metadata or summaries.

## Known Weak Points

1. **Formula constants embedded in scripts**

   The Bose critical-temperature prefactor `0.94` is a notebook formula constant. It is not currently inconsistent, but it is not visible in `configs/notebook_v1_defaults.json`.

2. **Numerical solver bracket embedded in scripts**

   The self-consistent total-atom root-solver bounds `[3e4, 1e7]` are procedural. They should remain harmless for the current notebook default, but they are not externally documented in config.

3. **Tests duplicate current outputs**

   Recovery tests protect against accidental drift, but some expected values are frozen from current config-derived outputs. If notebook defaults change intentionally, those tests must be refreshed deliberately.

4. **Exact RNG equivalence remains limited**

   Stochastic recovery is reproducible with explicit seeds, but exact historical notebook-global RNG call-order equivalence is not guaranteed. This is correctly documented as statistical/noise-model recovery rather than exact pixel-level notebook replay.

## Recommendations

No code or physics changes are required before continuing notebook-aligned recovery.

Recommended small follow-up, if traceability is tightened later:

1. Add a `formula_constants` or `numerical_procedure` section to `configs/notebook_v1_defaults.json` for:
   - Bose critical-temperature prefactor `0.94`
   - self-consistent total-atom solver bracket `[3e4, 1e7]`
2. Keep stage-specific values explicitly separated:
   - single-frame camera settings
   - multishot pulse settings
   - display-only figure settings
3. Keep stochastic recovery seed-based, and continue avoiding claims of exact notebook-global RNG replay unless call order is fully reconstructed.

## Conclusion

The notebook-aligned recovery pipeline is parameter-consistent for the current Version 1 defaults. The config file is functioning as the canonical source for active recovery parameters, metadata matches the configured values, and the remaining hard-coded values are either display/procedure constants or notebook formula coefficients rather than hidden alternate physics.

The main improvement opportunity is documentation/config traceability for formula and solver constants, not a correction to simulator behaviour.
