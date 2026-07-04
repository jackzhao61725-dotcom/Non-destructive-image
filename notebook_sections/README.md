# Notebook Section Exports

These files are a mechanical first pass at converting `../1 calculations revised 2  multishot  6  extended.ipynb` into source-controlled Python section files. They are intended to make the existing notebook easier to review and refactor without changing the scientific implementation.

Important constraints:

- The notebook remains the reference implementation.
- Equations, numerical algorithms, parameter values, and figure logic should remain unchanged.
- Code was copied from notebook cells into logical section files; notebook magics such as `%matplotlib inline` are commented out so the exports are valid Python text.
- Some cells appear in more than one section when they currently mix responsibilities, for example camera helpers embedded in imaging/demo cells. Later refactors should deduplicate these mechanically after regression checks.
- These section files are an intermediate refactoring aid, not yet a package API. Execute the original notebook for authoritative results until dependency order has been normalised.

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
