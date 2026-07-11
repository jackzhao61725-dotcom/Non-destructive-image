# Zenodo Release Checklist

This checklist prepares a future repository archive. It does not create a
release and does not mint a DOI by itself.

## Pre-Release Checks

1. Finalise the approved `main` branch.
2. Confirm these files exist:
   - `README.md`
   - `CITATION.cff`
   - `.zenodo.json`
   - `scripts/run_all_dissertation_figures.py`
   - `docs/reproducibility.md`
   - `docs/figure_index.md`
3. Choose and add a license before public archival or reuse.
4. Run validation:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
python scripts\run_all_dissertation_figures.py
```

5. Inspect:

```text
results/reproducibility_manifest.json
```

6. Confirm no unapproved generated outputs or feature-branch-only figures are
   included.

## GitHub Release

1. Create a semantic version tag, for example:

```text
v1.0.0
```

2. Create a GitHub release from that tag.
3. Include a concise release description:
   - Version 1 notebook-aligned simulator;
   - dissertation figure reproduction workflow;
   - representative uncalibrated outputs;
   - calibration limitations.

## Zenodo Archive

1. Connect the GitHub repository to Zenodo through the user's Zenodo account.
2. Enable archiving for the repository.
3. Archive the GitHub release on Zenodo.
4. Wait for Zenodo to mint the DOI.
5. Copy the DOI badge into `README.md`.
6. Add the DOI to `CITATION.cff` and README.
7. If needed, create a small follow-up GitHub release that includes the DOI
   badge and final citation metadata.

## Do Not Claim Early

Do not claim a DOI exists before Zenodo has minted it.

Do not claim calibrated experimental prediction unless experimental
absorption / RAI calibration has actually been applied and documented.

Do not change `kappa_F` from the Version 1 placeholder inside release metadata
unless a specific calibration milestone has been completed.
