# PCI Orchestration Milestone Report

## Objective

Migrate only the PCI-specific scalar imaging orchestration into a small helper
while preserving the current notebook-equivalent PCI behaviour.

## PCI Orchestration Migrated

Added:

```text
simulate_pci_image(...)
```

in:

```text
src/non_destructive_image/imaging.py
```

The helper performs only the PCI-specific orchestration above the existing
shared Fourier core:

```text
object_field = exp(1j * phase_map)
pci_reference_field = t_p * exp(1j * theta)
image_field = simulate_fourier_image(object_field, pupil, pci_reference_field, return_intensity=False)
pci_image_intensity = abs(image_field)**2
```

It can optionally return intermediate fields for regression tests:

- `object_field`
- `scattered_field`
- `propagated_scattered_field`
- `pci_reference_field`
- `pci_image_intensity`

## Notebook Section Status

`notebook_sections/04_pci.py` was left unchanged.

This was the lower-risk choice for this milestone because the section remains a
notebook export/migration surface with plotting, SNR, and downstream globals
interleaved. The helper is tested directly against the saved PCI/DGI imaging
baseline before any notebook section delegates to it.

## Behaviour Preservation

The helper preserves the notebook section 7.3 PCI conventions:

- scalar phase map convention,
- `object_field = exp(1j * phase_map)`,
- scattered field `object_field - 1`,
- FFT/pupil propagation through the existing `simulate_fourier_image(...)`,
- PCI reference field `t_p * exp(1j * theta)`,
- intensity calculation `abs(field)**2`.

No scalar phase, Thomas-Fermi, FFT, phase-plate, OD, or lower-layer helper API
was changed.

## Baseline Arrays Used

The regression test compares against:

```text
regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
```

Baseline arrays used:

- `scalar_phase_map_rad`
- `pupil`
- `object_field`
- `scattered_field`
- `propagated_scattered_field`
- `pci_reference_field`
- `pci_image_intensity`
- `metadata_json`

## Tests Added

Added:

```text
tests/regression/test_pci_orchestration.py
```

The test verifies:

- object field matches the baseline,
- scattered field matches the baseline,
- propagated scattered field matches the baseline,
- PCI reference field matches the baseline,
- PCI image intensity matches the baseline,
- returned intensity shape is correct,
- returned intensity is finite and non-negative.

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

DGI orchestration was not introduced.

Faraday orchestration was not introduced.

Camera, noise, shot-noise, SNR, optimisation, and multi-shot migration were not
touched.

No physics equations were changed.

No frozen lower-layer helper APIs were changed.

`simulate_fourier_image(...)` was not changed.

