# Notebook-Aligned Unit Consistency Audit

Audited source state: `main` at `60af9ac` (`Audit notebook-aligned parameter consistency`).

This audit checks the current notebook-aligned recovery pipeline for internal unit consistency and verifies the dissertation-facing parameter table supplied for the MSc report.

Generated audit files:

- `results/notebook_aligned_recovery/unit_inventory.csv`
- `results/notebook_aligned_recovery/thesis_parameter_table_check.csv`
- `results/notebook_aligned_recovery/thesis_parameter_table_check.json`

No simulator code, helper APIs, notebook sections, existing baselines, existing recovery outputs, Faraday optimisation results, or deliverables were changed.

## Overall Conclusion

No unit-conversion bug was found in the notebook-aligned recovery pipeline.

The dissertation parameter table values are consistent with recovered notebook-aligned values. All checked differences are explained by rounding in the dissertation-facing display values.

The main risks are documentation/display risks:

- Notebook-style display conventions include density in `cm^-3` and column density in `cm^-2`, while the supplied thesis parameter table uses `m^-3` and `m^-2`.
- Column-density subscripts `x`, `y`, and `z` mean integration along that axis, not a lineout plotted along that coordinate.
- A multishot figure label uses `frame number s`; in the notebook this `s` is a frame index symbol, not seconds.

These are label/interpretation risks, not physics or numerical mismatches.

## Internal Unit Convention

The recovery pipeline uses SI internally unless explicitly documented otherwise:

| Quantity | Internal unit | Display unit used where appropriate |
| --- | --- | --- |
| atomic mass | kg | amu only as source descriptor |
| scattering length | m | bohr radii in config/source descriptor |
| trap frequencies | Hz in config | Hz as `omega / 2pi` |
| angular trap frequencies | rad s^-1 | rad s^-1 only when explicitly needed |
| harmonic oscillator length | m | um |
| chemical potential | J | nK via `mu / k_B` |
| temperature-equivalent chemical potential | K | nK |
| Thomas-Fermi radii | m | um |
| 3D density | m^-3 | m^-3 for thesis table; cm^-3 for notebook-style display |
| column density | m^-2 | m^-2 for thesis table; cm^-2 for notebook-style display |
| grid coordinates | m | um |
| detuning | Hz | GHz for compact labels |
| natural linewidth | rad s^-1 | rad s^-1 |
| resonant cross section | m^2 | m^2 |
| scalar phase | rad | rad |
| Faraday rotation | rad | rad |
| optical power | mW in config | mW |
| intensity | W m^-2 | W m^-2 |
| exposure / pulse duration | s internally | us for labels |
| photon number | photons or detected electrons | photons/pixel or electrons |
| camera normalised image | dimensionless `I/I0` | dimensionless |
| read noise | electrons RMS | e- |
| SNR | dimensionless | dimensionless |
| atom number | atom count | atom count |
| loss fraction | fraction | fraction or percent |

## Frequency Convention

Trap frequencies in `configs/notebook_v1_defaults.json` are stored in Hz:

```text
[293.0, 14.0, 233.0] Hz
```

The Thomas-Fermi calculation converts them to angular frequencies with:

```text
omega_i = 2 * pi * f_i
```

The recovered geometric mean is:

```text
omega_bar / 2pi = 98.50324226832537 Hz
```

This matches the thesis display value `98.5 Hz` by rounding.

Detuning is stored as Hz for the scalar phase stages:

```text
Delta = 1.5e9 Hz
```

The natural linewidth is stored as rad s^-1:

```text
Gamma = 185353966.5617978 rad s^-1
```

The notebook dimensionless detuning convention is:

```text
delta = 2 * Delta_Hz * 2*pi / Gamma_rad_s
```

The recovered value is:

```text
delta = 101.69491525423727
```

This convention is implemented consistently in both the recovery scripts and `src/non_destructive_image/light_atom.py`.

## Length and Density Convention

Lengths are stored internally in metres. Figure axes and lineout CSV files convert coordinates to micrometres with:

```text
x_um = x_m * 1e6
```

The recovered TF radii are:

```text
R_x = 1.1857996475343457 um
R_y = 24.817092623397375 um
R_z = 1.4911557799466235 um
```

These match the thesis display values `(1.2, 24.8, 1.5) um` by rounding.

The recovered peak density is:

```text
n_0 = 3.400213390025543e20 m^-3
```

The recovery outputs also include notebook-style `cm^-3` display values using:

```text
cm^-3 = m^-3 * 1e-6
```

Column densities are stored internally in `m^-2`:

```text
n_tilde_x = 5.375962452578468e14 m^-2
n_tilde_y = 1.125112141861065e16 m^-2
n_tilde_z = 6.760330466117988e14 m^-2
```

Notebook-style `cm^-2` values use:

```text
cm^-2 = m^-2 * 1e-4
```

For the thesis parameter table, use the SI `m^-3` and `m^-2` values, not the notebook display `cm` values.

## Axis Convention

The config records:

```text
imaging_axis = 0
axis_labels = ["x", "y", "z"]
```

For Stage 18.1, the plotted image is the column-density map after integrating along `x`, so the displayed transverse plane is `y,z`.

The symbols `n_tilde_x`, `n_tilde_y`, and `n_tilde_z` should be interpreted as peak column-density scalars integrated along `x`, `y`, and `z` respectively. They are not lineouts along those axes and should not be used as labels for full 2D column-density distributions.

For full column-density maps, use distribution notation:

```text
x-integrated map in the y,z plane -> n_col(y,z)
y-integrated map in the x,z plane -> n_col(x,z)
z-integrated map in the x,y plane -> n_col(x,y)
```

This is the main axis-convention risk for dissertation text and captions.

## Thesis Parameter Table Check

| Quantity | Recovered dissertation value | Reference value | Status |
| --- | ---: | ---: | --- |
| `omega_bar / 2pi` | `98.50324226832537 Hz` | `98.5 Hz` | `ROUNDING_ONLY` |
| `a_ho` | `0.786220367322207 um` | `0.79 um` | `ROUNDING_ONLY` |
| `mu / k_B` | `47.573107349836604 nK` | `47.6 nK` | `ROUNDING_ONLY` |
| `n_0` | `3.400213390025543e20 m^-3` | `3.4e20 m^-3` | `ROUNDING_ONLY` |
| `R_x` | `1.1857996475343457 um` | `1.2 um` | `ROUNDING_ONLY` |
| `R_y` | `24.817092623397375 um` | `24.8 um` | `ROUNDING_ONLY` |
| `R_z` | `1.4911557799466235 um` | `1.5 um` | `ROUNDING_ONLY` |
| `n_tilde_x` | `5.375962452578468e14 m^-2` | `5.4e14 m^-2` | `ROUNDING_ONLY` |
| `n_tilde_y` | `1.125112141861065e16 m^-2` | `1.1e16 m^-2` | `ROUNDING_ONLY` |
| `n_tilde_z` | `6.760330466117988e14 m^-2` | `6.8e14 m^-2` | `ROUNDING_ONLY` |

The full machine-readable check is in:

- `results/notebook_aligned_recovery/thesis_parameter_table_check.csv`
- `results/notebook_aligned_recovery/thesis_parameter_table_check.json`

## Hidden Conversion Factors Reviewed

The following conversion factors are present and were checked:

- `2*pi`: Hz to angular frequency for trap frequencies.
- `2 * Delta_Hz * 2*pi / Gamma_rad_s`: notebook dimensionless detuning convention.
- `2*pi*hbar`: reconstructs Planck's constant for photon energy.
- `1e6`: metres to micrometres.
- `1e9`: kelvin to nanokelvin and Hz to GHz display.
- `1e-6`: `m^-3` to `cm^-3`.
- `1e-4`: `m^-2` to `cm^-2`.
- `1e-3`: mW to W.
- `1e-6`: microseconds to seconds in multishot recovery.

No incorrect conversion direction was found.

## Stage-Specific but Valid Conventions

Some stage-specific settings are intentional:

- Single-frame deterministic/noisy camera recovery uses `P = 2.0 mW` and `tau = 100 us`.
- Multishot recovery uses `P = 3.5 mW` and `tau = 40 us`.
- Filmstrip photon counts use the multishot operating point, so `photons_per_pixel = 1072.6265629146646`, not the single-frame camera value `1532.3236613066642`.
- DGI uses `OD = 4.0` as a DGI-only reference parameter.
- Faraday uses `kappa_F = 1.0` as the Version 1 phenomenological placeholder.

These are not unit inconsistencies.

## Risks and Classifications

| Risk | Classification | Recommendation |
| --- | --- | --- |
| Notebook display density units are `cm^-3` / `cm^-2`, while thesis table uses `m^-3` / `m^-2` | `UNIT_LABEL_RISK` | In dissertation tables, explicitly label density as `m^-3` and column density as `m^-2`. |
| Column-density subscript may be confused with plotted lineout axis or 2D map label | `AXIS_CONVENTION_RISK` | Use `n_col(y,z)` for an x-integrated distribution; reserve `n_tilde_x` for its peak scalar value. |
| Multishot label `frame number s` could be read as seconds | `UNIT_LABEL_RISK` | In dissertation figures, use `frame index` or `shot index` instead. |
| Bose critical-temperature prefactor `0.94` and solver bounds remain script constants | `NEEDS_CONFIG_TRACEABILITY` | Same recommendation as the parameter audit: optionally move formula/procedure constants into config later. |

No `POSSIBLE_BUG`, `UNIT_MISMATCH`, or value mismatch was found.

## Recommended Next Steps

No code fix is required.

Before using these quantities in the dissertation:

1. Use the SI density and column-density values from `thesis_parameter_table_check.csv`.
2. In captions, state that Stage 18.1 is an `x`-integrated column-density map shown in the `y,z` plane.
3. Avoid using `cm^-3` or `cm^-2` table values unless the table explicitly asks for cgs units.
4. Rename `frame number s` to `frame index` in any future dissertation-ready plotting script.

## Plot Label Unit Conventions

Notebook-aligned recovery figures should use Matplotlib mathtext-compatible
labels rather than plain-text superscripts or copied word-processor equation
artifacts.

Use:

- `$\mu\mathrm{m}$` for micrometres in plotted axes.
- `$\mathrm{m}^{-3}$` and `$\mathrm{m}^{-2}$` for SI density and column-density units.
- `$\mathrm{cm}^{-3}$` and `$\mathrm{cm}^{-2}$` only for notebook-style display quantities that are explicitly converted from SI.
- `$n_{\mathrm{col}}(a,b)$` for a full 2D column-density distribution in the displayed `a,b` plane.
- `$\tilde{n}_i$` only for the peak column-density scalar integrated along axis `i`.
- `$\bar{\omega}/2\pi$` for displayed trap frequency in Hz.
- `$\mu/k_B$` for chemical-potential temperature equivalent.
- `$\phi$ (rad)` for scalar phase.
- `$\theta_F$ (rad)` for Faraday rotation.
- `frame index` or `shot index` for multishot frame counters.

Avoid:

- `m^-2`, `m^-3`, `m^(-2)`, or `m^(-3)` in figure labels.
- mixed `um`, Unicode `μm`, and copied text artifacts in the same figure family.
- using `s` in axis labels where it means a frame index rather than seconds.
- unlabeled density conversions between SI and notebook-style cgs display units.

The notebook-aligned recovery scripts now share common labels through
`scripts/plot_label_utils.py` so regenerated SVGs use consistent unit rendering.
