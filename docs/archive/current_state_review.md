# Current Code State Review

This review summarises the repository after the non-destructive notebook audit, section export, helper extraction, and validation work. It is intentionally conservative: the original notebook is still the scientific reference implementation.

## 1. Reference status

- Reference notebook: `1 calculations revised 2  multishot  6  extended.ipynb`.
- The notebook itself has not been edited during this refactor sequence.
- The exported files and helper package are intended to support future cleanup, not to replace the notebook yet.

## 2. What currently works

### Notebook exports

- `notebook_sections/` contains a mechanical export of the notebook into logical 00-10 section files.
- `scripts/export_notebook_sections.py` can regenerate those files from the notebook JSON.
- `scripts/validate_notebook_sections.py` confirms the checked-in section files are exactly in sync with the generated output and are Python-syntax valid.

Current validation result in this environment:

```text
Notebook section exports are in sync and syntactically valid.
```

### Helper package

The helper package under `src/non_destructive_image/` now contains small, explicit-parameter extractions of repeated notebook patterns:

- `atomic_model.py` — Thomas-Fermi state and recoil formulas.
- `profiles.py` — 2D Thomas-Fermi column-profile expression.
- `fourier.py` — FFT/pupil propagation helper.
- `camera.py` — camera binning, shot/read-noise, and normalisation helpers.
- `light_atom.py` — detuning, scalar phase, residual optical depth, intensity, and scattering formulas.

These helpers are deliberately not wired back into the notebook yet. They are ready to be used in later small refactor PRs, one replacement at a time.

### Tests and checks

- Python syntax compilation succeeded for the helper package, tests, and scripts.
- Export validation succeeded.
- The pytest suite is present, but in this execution environment it skipped because `numpy` is not installed.
- `requirements.txt` documents the intended runtime/test dependencies: `numpy` and `pytest`.

## 3. Equivalence to the old notebook

The current repository is equivalent to the old notebook in the following limited sense:

- The notebook remains unchanged and therefore still represents the old working implementation.
- The section files are mechanical text exports from the notebook and validate as in sync with the notebook JSON.
- The helper functions are direct algebraic extractions of repeated notebook expressions and are covered by equivalence tests that compare against those expressions.

The current repository is not yet a full replacement for executing the notebook because:

- The section exports preserve notebook cells, including mixed responsibilities and dependencies.
- The helper package has not yet been wired into the exported section files or the notebook.
- Full numerical/figure regression against a notebook execution was not run in this environment because the scientific Python dependency stack is unavailable here.

## 4. Current file map

| Path | Purpose | Status |
|---|---|---|
| `1 calculations revised 2  multishot  6  extended.ipynb` | Original scientific notebook | Reference implementation; unchanged |
| `README.md` | Project overview and validation commands | Current |
| `docs/notebook_audit.md` | Human-readable audit and refactor plan | Current |
| `docs/notebook_audit.json` | Machine-readable cell/function inventory | Current |
| `docs/current_state_review.md` | This state review | Current |
| `notebook_sections/` | Mechanical logical section exports | Generated and validated |
| `scripts/export_notebook_sections.py` | Regenerates section exports | Current |
| `scripts/validate_notebook_sections.py` | Checks exports are in sync and compile | Current |
| `src/non_destructive_image/` | Extracted helper package | Initial helpers present; not yet wired into notebook |
| `tests/test_helpers.py` | Helper equivalence tests | Ready; requires `numpy` to execute fully |
| `requirements.txt` | Minimal Python dependencies | Current |
| `pytest.ini` | Pytest configuration | Current |
| `.gitignore` | Python/cache ignore rules | Current |

## 5. Recommended next steps

1. Install the dependency stack in an environment with package access: `numpy` and `pytest` at minimum.
2. Run `python3 scripts/validate_notebook_sections.py`.
3. Run `python3 -m pytest -q`.
4. Execute the original notebook once to capture baseline scalar outputs and representative arrays.
5. Wire helpers into exported section files one group at a time, starting with camera binning/noise.
6. After each replacement, compare against the captured baseline before moving to the next group.

## 6. Bottom line

We have enough structure to begin careful behaviour-preserving refactoring, but the old notebook should still be treated as the only authoritative runnable scientific workflow until the helper package is wired in and full numerical regression has been run in a complete scientific Python environment.
