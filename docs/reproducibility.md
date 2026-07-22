# Reproducibility

- **Status:** active environment, command and output-retention authority
- **Active consumers:** local verification, maintained figure generation and
  sealed-result inspection
- **Update trigger:** a change to the supported interpreter, dependencies,
  configs, generators, default run plan, retained outputs or verification gates
- **Retirement rule:** replace only with one indexed successor covering both
  current and historical reproduction boundaries

## Environment

This worktree has a repository-local Python 3.12 `.venv` with the project
installed in editable mode. From the repository root, use its explicit
interpreter so Windows aliases and other editable clones cannot affect a run:

```powershell
$repoRoot = (Resolve-Path -LiteralPath ".").Path
$projectPython = "$repoRoot\.venv\Scripts\python.exe"
$env:PYTHONPATH="$repoRoot\src;$repoRoot"
$env:PYTHONUTF8="1"
& $projectPython -m pip check
```

If the environment must be rebuilt, create it with a verified Python 3.12
interpreter, then run `& $projectPython -m pip install -e '.[dev]'`. Do not use
the current WindowsApps `python` alias or an interpreter without the project
dependencies.

Before any numerical run, verify that the selected interpreter imports this
checkout rather than another editable clone:

```powershell
& $projectPython -c "import sys, non_destructive_image; print(sys.executable); print(non_destructive_image.__file__)"
```

The package path must lie under the repository from which the command is being
run. On Windows, do not assume that a bare `python` or `pytest` resolves to the
intended environment; use the environment's explicit interpreter and
`python -m ...` form.

## Validation

Run the full test suite:

```powershell
& $projectPython -m pytest -q
```

## One-command regeneration

Inspect the ordered plan:

```powershell
& $projectPython scripts\run_all_dissertation_figures.py --dry-run
```

Regenerate the maintained default outputs:

```powershell
& $projectPython scripts\run_all_dissertation_figures.py
```

The script writes `results/reproducibility_manifest.json`, containing the
Python environment, Git commit, commands, configs and expected outputs.
The default run excludes frozen Figures 5.2 and 5.4. They are listed as blocked
pending items and their canonical directories cannot be overwritten by the
Version 1 generators. Allow several minutes for the retained notebook-aligned
and current single-frame outputs.

## Active configs

| Role | Config |
| --- | --- |
| historical notebook regression | `configs/notebook_v1_defaults.json` |
| current optical/detector and implemented Version 1 sequence reference | `configs/dissertation_v3_orca_fusion.json` |
| shared detector-dependent scans | `configs/dissertation_plots_v2_orca_fusion.json` |
| figures | `configs/figure_4_2.json`, `figure_5_1.json`, `figure_5_2.json`, `figure_5_4.json` |
| active morphology inverse | `configs/reconstruction_morphology_benchmark_v4_orca_fusion_m10.json` |
| inverse credibility study | `configs/reconstruction_credibility_v2_orca_fusion_m10.json` |
| held-out physical observables | `configs/reconstruction_observables_v1_orca_fusion_m10.json` |
| planned seven-condition historical-inverse replay; no accepted result | `configs/dpfi_initial_condition_suite_v1_orca_fusion_m10.json` |

## Dependent outputs

The frozen Version 1 Figure 5.4 uses the same physical and image-quality
contract as the frozen Version 1 Figure 5.2.
Changing magnification, QE, read noise, numerical aperture, the atomic Faraday
conversion, the effective apparatus response or the disturbance model requires
regeneration of the dependent Chapter 5 outputs and their regression tests.

## Focused current commands

```powershell
& $projectPython scripts\generate_detuning_tradeoff_plot.py --config configs\dissertation_plots_v2_orca_fusion.json
& $projectPython scripts\generate_figure_4_2.py --config configs\figure_4_2.json
& $projectPython scripts\generate_figure_5_1.py --config configs\figure_5_1.json
```

Figures 5.2 and 5.4 are retained only as frozen Version 1 screening outputs.
Their generators require an explicit non-canonical `--output-dir` for a
historical reproducibility check and refuse to overwrite the retained result
directories. Their next dissertation-facing regeneration must use the approved
shared heating replacement.

## Reconstruction workflow

The active benchmark and credibility scripts default to the ORCA-Fusion/M10
configs:

```powershell
& $projectPython scripts\generate_reconstruction_morphology_benchmark.py
& $projectPython scripts\plot_reconstruction_morphology_benchmark.py
& $projectPython scripts\run_reconstruction_curvature_range_check.py --config configs\reconstruction_curvature_range_check_v2_orca_fusion_m10_dual_port.json
& $projectPython scripts\run_reconstruction_curvature_range_check.py --config configs\reconstruction_curvature_range_check_v2_orca_fusion_m10_dark_field.json
& $projectPython scripts\run_reconstruction_credibility_study.py
& $projectPython scripts\generate_reconstruction_observable_benchmark.py --validate-only
& $projectPython scripts\generate_reconstruction_observable_benchmark.py
& $projectPython scripts\plot_reconstruction_observable_benchmark.py
```

The curvature command deliberately requires an explicit config because the two
readouts have separate range-check records. The morphology search and
credibility ensemble are expensive and are kept outside
`run_all_dissertation_figures.py`. The observable generator first verifies the
sealed morphology run, then replays only its 60 frozen held-out fits. It took
1715.90 seconds on the recorded runtime. Use `--validate-only` to check the
source run, trial axes and object-space support without fitting; the plot command
only verifies and renders the sealed observable tables. Existing results should
be inspected and verified before regeneration.

## Output retention

The figure and data index is the authority for which result families have a
current consumer. A maintained output needs an explicit config, a deterministic
generator, machine-readable provenance and a regression check when it supports
a dissertation claim. Trial renders, alternate crops and one-off scans do not
remain in the repository after review. Recover a deleted artifact from Git
history or regenerate it from recorded inputs; do not create an archive
directory.

## Provenance requirements

Every retained output must record:

- source config paths;
- detuning, power, exposure and imaging axis;
- optical sampling and detector parameters;
- signal normalisation and estimator/ROI;
- sequence and stopping model;
- random seed when a stochastic frame is rendered;
- calibration status;
- output path and Git commit.

Generated camera images are stochastic illustrations. Reported SNR is computed
from expected counts and the analytic noise variance.

## Scientific boundary

The sealed reconstruction result set is internally reproducible under its
frozen `kappa_F=1` operator, but that conversion has been superseded. It is
method-development evidence, not a current `166Er` amplitude, SNR or sequence
prediction. The stored Figures 5.2 and 5.4 use signed `kappa_F=-45/91` but retain
the frozen Version 1 heating model. New quantitative sequence outputs must also
use the approved heating replacement and still require an effective apparatus
calibration. Validation requires separate data not used to set those parameters.
