# Repository Large-File Hygiene Review

## Current File Sizes

Measured after updating the bundle policy:

```text
3.76 MB   1 calculations revised 2  multishot  6  extended.ipynb
2.86 MB   deliverables/non_destructive_image_code_bundle.zip
21.89 MB  regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
76.43 MB  regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

Files above 10 MB at the current repository tip:

```text
76.43 MB  regression/baseline/imaging/faraday_imaging_baseline_v1.npz
21.89 MB  regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
```

Before this hygiene update, the generated deliverable bundle was about
98.64 MB because it embedded the large numerical `.npz` baselines.

## Risk Assessment

GitHub warns at 50 MB and rejects single files above 100 MB.

The Faraday baseline is the highest current risk. At about 76.43 MB it is below
GitHub's 100 MB hard limit, but close enough that future additions to the same
artifact format could exceed the limit.

The PCI/DGI baseline is about 21.89 MB. It is not close to the hard limit, but
it contributes to repository weight.

The original notebook is about 3.76 MB and is not a large-file risk.

The regenerated deliverable bundle is now about 2.86 MB and is no longer close
to GitHub's warning or hard-limit thresholds.

Because this fix does not rewrite Git history, older commits still contain the
previous larger bundle object. This is intentional for a conservative no-history-
rewrite recovery path. The current repository tip no longer generates or tracks
a large bundle artifact.

## Chosen Policy

The generated code bundle is a portable source archive, not a container for
large numerical baselines.

`scripts/create_code_bundle.py` now excludes generated `.npz` files under:

```text
regression/baseline/
```

The repository continues to track the current `.npz` baselines so existing
regression tests keep working.

The bundle README explains that large baselines are kept separately in the
source repository and can be regenerated with:

```bash
python scripts/generate_pci_dgi_imaging_baseline.py
python scripts/generate_faraday_imaging_baseline.py
```

## Future Recommendations

If baseline files continue to grow, consider Git LFS for generated numerical
artifacts.

If future baselines only need representative regression coverage, consider
using smaller deterministic grids for baseline artifacts. For the Faraday
baseline, a smaller grid would likely reduce file size substantially, but that
would change the stored regression artifact and should be a separate explicit
milestone.

For releases, prefer attaching generated bundles and large baselines as GitHub
release artifacts rather than repeatedly committing regenerated archives.

## Validation

The hygiene update changes only packaging policy and documentation.

No physics equations were changed.

No helper APIs were changed.

No notebook logic was changed.

No simulator behavior was changed.
