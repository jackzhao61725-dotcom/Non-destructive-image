# Scientific regression baselines

- **Status:** maintained historical regression-data contract
- **Active consumers:** notebook-output and deterministic imaging regression tests
- **Update trigger:** an intentional change to a recorded notebook convention,
  baseline generator, metadata schema or comparison tolerance
- **Retirement rule:** remove a family only with its generator, tests and all
  provenance references in the same change

This directory freezes selected outputs and numerical conventions from the
historical reference notebook. The files test whether maintained code remains
compatible with a recorded convention. A passing regression does not establish
that the historical convention is the current physical authority.

## Maintained baseline families

- `notebook_outputs.json` records text and hashes of rich outputs already stored
  in `1 calculations revised 2  multishot  6  extended.ipynb`. It is extracted
  without executing the notebook by
  `scripts/extract_notebook_output_baseline.py`.
- `imaging/pci_dgi_imaging_baseline_v1.npz` contains deterministic PCI and DGI
  arrays produced by `scripts/generate_pci_dgi_imaging_baseline.py`. The
  generator mirrors the notebook Section 7.3 / cell 18 convention; it does not
  execute the full notebook.
- `imaging/faraday_imaging_baseline_v1.npz` contains deterministic Faraday and
  camera arrays produced by `scripts/generate_faraday_imaging_baseline.py`. The
  generator mirrors the historical Section 17.2 / cell 51 convention; it does
  not execute the full notebook. It preserves the historical `kappa_F = 1`
  placeholder for compatibility and must not be read as the current physical
  coefficient for 166Er.

## File and metadata contract

Numerical families use compressed NumPy archives created with
`numpy.savez_compressed`. Arrays have descriptive names with units where
needed, and structured metadata is stored as the `metadata_json` string array.
Archives must not contain pickled Python objects. Each numerical generator
records the source notebook path, its SHA-256 hash, the relevant notebook
location, whether the notebook was executed, and the Python and NumPy versions.

The current imaging baselines are deterministic and contain no sampled camera
noise. Any future stochastic baseline must use a fixed, documented random seed
and generator type. Store both in its metadata. Compare deterministic pre-noise
arrays across environments; do not compare independent stochastic draws.

Regression tests compare floating-point imaging arrays with
`rtol=1e-10, atol=1e-12`. Algebraic scalar extractions may use
`rtol=1e-12, atol=0` when the implementation and backend make that meaningful.
Any relaxation must be justified in the test that requires it. Figures are not
regression targets; compare the arrays from which they are drawn.

## Update rule

Do not edit baseline data by hand or regenerate it merely to make a failing test
pass. First determine whether the failure is an implementation regression or an
intentional change to the recorded notebook convention. An intentional update
must preserve an origin independent of the code path being tested: extract the
stored notebook outputs or update the explicit notebook-equivalent generator,
then regenerate through the corresponding maintained script:

Use the explicit `$projectPython` defined in
[the reproducibility guide](../../docs/reproducibility.md):

```powershell
& $projectPython scripts\extract_notebook_output_baseline.py
& $projectPython scripts\generate_pci_dgi_imaging_baseline.py
& $projectPython scripts\generate_faraday_imaging_baseline.py
```

Review the source-notebook hash, metadata, array changes, and affected
regression tests together. Add a new baseline family only when it has both an
implemented generator and an active test or analysis consumer.
