# Colab Validation Guide

This guide explains how to validate the repository in Google Colab or another complete scientific Python environment. The current lightweight execution environment may skip NumPy-backed tests when NumPy is unavailable; Colab should run those tests normally once dependencies are installed.

## Validation status to confirm

For the shared PCI/DGI Fourier-imaging core:

> Scope accepted; architecture accepted; physics boundary accepted; numerical acceptance pending full scientific-environment test execution.

The milestone should only be fully closed after the NumPy-backed tests pass in Colab or another scientific Python environment.

## 1. Clone the repository

Replace `<repo-url>` with the GitHub URL for this repository.

```bash
git clone <repo-url>
cd <repo>
```

If validating a pull-request branch, check out that branch before installing dependencies.

## 2. Install dependencies

This repository currently uses `requirements.txt` for the minimal test dependencies. Install it first, then install the extra scientific/notebook packages needed for full notebook-style validation.

```bash
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install scipy matplotlib
```

The project does not currently define a packaging file such as `pyproject.toml` or `setup.py`, so use `PYTHONPATH=src` when running tests and scripts instead of `pip install -e .`.

## 3. Run the validation suite

```bash
PYTHONPATH=src pytest -q
python scripts/validate_notebook_sections.py
```

Alternatively, run the convenience script:

```bash
python scripts/run_validation.py
```

The convenience script sets `PYTHONPATH=src` for pytest and then runs the notebook-section validation script.

## 4. Expected result

- `pytest -q` should run the NumPy-backed helper tests rather than skipping them.
- The PCI/DGI Fourier-core test should compare `simulate_fourier_image(...)` against the original notebook-equivalent FFT expression.
- `python scripts/validate_notebook_sections.py` should report that unmigrated exports are in sync and migrated sections are syntactically valid.

## 5. If tests skip unexpectedly

If NumPy-backed tests still skip in Colab, verify:

```bash
python - <<'PY'
import numpy
print(numpy.__version__)
PY
```

If that import fails, reinstall dependencies with:

```bash
python -m pip install -r requirements.txt scipy matplotlib
```

Do not treat the Imaging Milestone 1 numerical acceptance as complete until the NumPy-backed tests run and pass.
