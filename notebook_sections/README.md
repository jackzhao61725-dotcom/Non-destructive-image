# Notebook Section Exports

- **Status:** historical mechanical exports retained for regression and provenance
- **Active consumers:** the deterministic export validator, notebook-recovery
  scripts and notebook-aligned regression tests
- **Update trigger:** the source notebook changes or the deterministic export
  validator changes its accepted cell mapping
- **Retirement rule:** remove only after every active regression and provenance
  consumer has migrated to another sealed source

These files are a mechanical export of
`../1 calculations revised 2  multishot  6  extended.ipynb`. They make the
historical notebook conventions reviewable without changing them. The
maintained package and active configs are the current implementation authority.

Important constraints:

- The notebook is a historical regression reference, not the active physical or
  software implementation.
- Equations, numerical algorithms, parameter values, and figure logic should remain unchanged.
- Code was copied from notebook cells into logical section files; notebook magics such as `%matplotlib inline` are commented out so the exports are valid Python text.
- Some cells appear in more than one section when they currently mix responsibilities, for example camera helpers embedded in imaging/demo cells. Later refactors should deduplicate these mechanically after regression checks.
- These section files are not a package API. Use the maintained recovery scripts
  and validators for historical reproduction; do not execute the notebook as a
  source of current dissertation results.

Section files:

- `00_imports.py` — 00 Imports; notebook cells [0, 1, 2].
- `01_parameters.py` — 01 Parameters; notebook cells [3, 4, 5, 6].
- `02_atomic_model.py` — 02 Atomic Model; notebook cells [7, 8].
- `03_light_atom_interaction.py` — 03 Light-Atom Interaction; notebook cells [9, 10, 11, 12, 13, 14, 32, 49].
- `04_pci.py` — 04 PCI; notebook cells [15, 16, 17, 18, 21, 22, 35, 36].
- `05_dgi.py` — 05 DGI; notebook cells [19, 20].
- `06_faraday.py` — 06 Faraday; notebook cells [48, 49, 50, 51, 52, 53, 54, 55, 84, 85, 86, 87, 88, 89, 90, 91].
- `07_camera.py` — 07 Camera; notebook cells [20, 44].
- `08_shot_noise.py` — 08 Shot Noise; notebook cells [21, 22, 23, 24, 25, 26, 33, 34].
- `09_multishot_simulation.py` — 09 Multi-shot Simulation; notebook cells [39, 40, 41, 42, 43, 44, 45, 46, 92, 93].
- `10_analysis.py` — 10 Analysis; notebook cells [27, 28, 29, 30, 31, 37, 38, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 94, 95].

Unassigned notebook cells: [47].
