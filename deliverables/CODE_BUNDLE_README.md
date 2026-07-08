# Non-destructive Image Code Bundle

This directory contains the generated code bundle for the current repository
state.

## Contents

`non_destructive_image_code_bundle.zip` is a portable archive of the project
source needed to inspect, validate, and continue the notebook refactor. It
contains:

- `README.md`
- `requirements.txt`
- `pytest.ini`
- the original notebook, `1 calculations revised 2  multishot  6  extended.ipynb`
- `docs/`
- `notebook_sections/`
- `src/`
- `scripts/`
- `tests/`
- `regression/`, excluding generated `.npz` numerical baseline artifacts

The archive excludes local development state such as `.git/`, `.venv/`,
Python caches, pytest caches, editor metadata, temporary files, and previously
generated zip files inside `deliverables/`.

Large generated numerical baselines are intentionally kept out of the portable
code bundle to keep the archive comfortably below GitHub's large-file warning
thresholds. The source repository tracks the current baseline `.npz` files
separately under `regression/baseline/`, and the bundle includes the scripts and
documentation needed to regenerate them.

## Scientific Source Of Truth

The original notebook remains the authoritative scientific implementation for
physics equations, numerical behavior, parameters, and figures.

The `src/non_destructive_image/` helper package is a refactored support layer
that extracts repeated notebook formulas and utility operations into small,
testable functions. It should be treated as support code for conservative
migration work, not as a replacement for the reference notebook unless future
validated migration steps explicitly say so.

## Install Dependencies

From the repository root, install the lightweight validation dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run Tests

On Windows PowerShell:

```powershell
$env:PYTHONPATH = "src;."
pytest -q
```

On Linux or macOS shells:

```bash
PYTHONPATH=src:. pytest -q
```

## Validate Notebook Section Exports

On Windows PowerShell:

```powershell
python scripts\validate_notebook_sections.py
```

On Linux or macOS shells:

```bash
python scripts/validate_notebook_sections.py
```

This validation checks that the source-controlled section exports match the
reference notebook export process and compile as Python.

## Regenerate Numerical Baselines

The source repository currently tracks large representative `.npz` baselines
for PCI/DGI and Faraday imaging regression tests. They are omitted from the code
bundle by default, but can be regenerated from the repository root with:

```bash
python scripts/generate_pci_dgi_imaging_baseline.py
python scripts/generate_faraday_imaging_baseline.py
```

After regenerating these files, the full regression suite can be run with
`pytest -q`.

## Regenerate The Bundle

If `scripts/create_code_bundle.py` is present, regenerate the archive from the
repository root with:

```bash
python scripts/create_code_bundle.py
```

The script writes:

```text
deliverables/non_destructive_image_code_bundle.zip
```

