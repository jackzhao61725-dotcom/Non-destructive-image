# PCI/DGI Imaging Baseline Milestone Report

## Objective

Generate deterministic PCI/DGI scalar imaging baseline arrays before introducing
PCI-specific or DGI-specific orchestration helpers.

This milestone is baseline generation only. It does not migrate notebook logic,
change physics, add optimisation, or introduce `simulate_pci_image(...)` or
`simulate_dgi_image(...)`.

## Arrays Saved

The baseline file is:

```text
regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
```

It contains:

- `grid_axis_m`
- `grid_a_m`
- `grid_b_m`
- `spatial_frequency_axis_m_inv`
- `spatial_frequency_x_m_inv`
- `spatial_frequency_y_m_inv`
- `radii_m`
- `column_density_m2`
- `phase_peak_rad`
- `thomas_fermi_profile`
- `scalar_phase_map_rad`
- `object_field`
- `pupil`
- `scattered_field`
- `propagated_scattered_field`
- `pci_reference_field`
- `pci_image_intensity`
- `dgi_reference_field`
- `dgi_image_intensity`
- `metadata_json`

## Parameters Used

The baseline uses the notebook section 7.3 reference PCI/DGI scalar imaging
setup:

- imaging axis: `x`
- transverse plane: `y,z`
- detuning: `1.5e9 Hz`
- grid size: `1024 x 1024`
- field of view: `100e-6 m`
- wavelength: `401e-9 m`
- probe diameter: `24e-3 m`
- first imaging lens focal length: `150e-3 m`
- numerical aperture: `(D_probe / 2) / f1`
- PCI plate amplitude transmittance: `t_p = 0.95`
- PCI phase shift: `theta = pi/2`
- DGI stop optical density: `OD = 4.0`

Atomic and optical parameters match the current reference notebook values for
the `166Er` 401 nm transition and BEC operating point.

## Generation Method

The generator is:

```text
scripts/generate_pci_dgi_imaging_baseline.py
```

It directly mirrors the current notebook-equivalent scalar PCI/DGI Fourier
imaging equations from notebook section 7.3 / cell 18:

```text
object_field = exp(1j * scalar_phase_map)
scattered_field = object_field - 1
propagated_scattered_field = ifft2(fft2(scattered_field) * pupil)
PCI field = t_p * exp(1j * theta) + propagated_scattered_field
DGI field = 10**(-OD/2) + propagated_scattered_field
intensity = abs(field)**2
```

This is a representative deterministic notebook-equivalent array baseline, not
a full notebook execution dump. No stochastic camera noise is included.

## Regeneration

From the repository root:

```bash
python scripts/generate_pci_dgi_imaging_baseline.py
```

This rewrites:

```text
regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
```

## Tests Added

The regression test is:

```text
tests/regression/test_pci_dgi_imaging_baseline.py
```

It checks:

- the baseline file exists,
- required keys are present,
- array shapes are as expected,
- arrays contain finite values,
- PCI and DGI intensity arrays are non-negative,
- metadata is present,
- regenerating the baseline produces arrays equal to the stored `.npz` values.

## Validation Results

Validation was run with:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
python scripts\create_code_bundle.py
```

Results:

```text
pytest -q: passed
notebook section validation: passed
code bundle generation: passed and remained reproducible
```

## Scope Confirmation

No physics equations were changed.

No helper APIs were changed.

No notebook logic was changed.

No simulator behaviour was changed.

No PCI-specific, DGI-specific, Faraday, camera, noise, shot-noise, SNR, or
multi-shot orchestration helper was introduced.

