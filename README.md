# Non-destructive image notebook refactor

This repository currently treats `1 calculations revised 2  multishot  6  extended.ipynb` as the scientific reference implementation for the ultracold atom imaging calculations.

The refactor is proceeding conservatively:

1. audit the notebook,
2. export notebook cells into source-controlled section files,
3. add validation checks so future refactors can be behaviour-preserving,
4. only then begin renaming variables and extracting shared helpers.

## Where to look

- `1 calculations revised 2  multishot  6  extended.ipynb` — reference notebook; authoritative for physics and figures.
- `docs/notebook_audit.md` — human-readable audit and refactoring plan.
- `docs/notebook_audit.json` — machine-readable notebook inventory.
- `notebook_sections/` — mechanical Python exports grouped into the planned logical sections.
- `scripts/export_notebook_sections.py` — regenerates the section exports from the notebook.
- `scripts/validate_notebook_sections.py` — checks that the exports are in sync and Python-syntax valid.

## Extracted helpers

The first behaviour-preserving helper extractions live in `src/non_destructive_image/`:

- `atomic_model.py` — Thomas-Fermi condensate-state and recoil helper formulas.
- `profiles.py` — Thomas-Fermi column-profile expression.
- `fourier.py` — FFT/pupil propagation pattern.
- `camera.py` — camera binning, noise, and normalisation helpers.
- `light_atom.py` — detuning, scalar phase, residual OD, intensity, and scattering helpers.

These helpers are not yet wired back into the notebook. They are covered by small equivalence tests in `tests/test_helpers.py` so later notebook edits can replace duplicated code safely.

## Validation

Install the notebook-refactor test dependencies listed in `requirements.txt` when network/package access is available. Then run the lightweight checks with:

```bash
python3 scripts/validate_notebook_sections.py
PYTHONPATH=src:. python3 -m pytest -q
```

The export validation does not execute the physics simulation. It verifies that the section files are exactly reproducible from the notebook and that the exported Python files compile. The pytest command checks the extracted helpers against the exact notebook expressions they replace.
