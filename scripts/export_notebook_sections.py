"""Export the research notebook into logical Python section files.

This script is intentionally mechanical: it copies notebook cells into section
files without changing equations, numerical algorithms, parameter values, or
plotting logic. Markdown cells are emitted as Python comments and IPython magics
are commented out so the generated files are valid Python text.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

NOTEBOOK_PATH = Path("1 calculations revised 2  multishot  6  extended.ipynb")
OUTPUT_DIR = Path("notebook_sections")

SECTIONS: list[tuple[str, str, Iterable[int]]] = [
    ("00_imports.py", "00 Imports", range(0, 3)),
    ("01_parameters.py", "01 Parameters", range(3, 7)),
    ("02_atomic_model.py", "02 Atomic Model", range(7, 9)),
    ("03_light_atom_interaction.py", "03 Light-Atom Interaction", [*range(9, 15), 32, 49]),
    ("04_pci.py", "04 PCI", [15, 16, 17, 18, 21, 22, 35, 36]),
    ("05_dgi.py", "05 DGI", [19, 20]),
    ("06_faraday.py", "06 Faraday", [*range(48, 56), *range(84, 92)]),
    ("07_camera.py", "07 Camera", [20, 44]),
    ("08_shot_noise.py", "08 Shot Noise", [21, 22, 23, 24, 25, 26, 33, 34]),
    ("09_multishot_simulation.py", "09 Multi-shot Simulation", [39, 40, 41, 42, 43, 44, 45, 46, 92, 93]),
    ("10_analysis.py", "10 Analysis", [27, 28, 29, 30, 31, 37, 38, *range(56, 84), 94, 95]),
]


def _comment_markdown(source: str) -> list[str]:
    return ["# " + line if line else "#" for line in source.splitlines()]


def _comment_ipython_magics(source: str) -> list[str]:
    lines: list[str] = []
    for line in source.splitlines():
        if line.lstrip().startswith("%"):
            lines.append("# " + line)
        else:
            lines.append(line)
    return lines


def _render_cell(cell_index: int, cell: dict) -> list[str]:
    source = "".join(cell.get("source", [])).rstrip("\n")
    lines = ["", f"# %% [cell {cell_index}: {cell['cell_type']}]"]
    if cell["cell_type"] == "markdown":
        lines.extend(_comment_markdown(source))
    else:
        lines.extend(_comment_ipython_magics(source))
    return lines


def _render_section(filename: str, title: str, cells: Iterable[int], notebook: dict) -> str:
    lines = [
        f"# {title}",
        "#",
        f"# Exported from: {NOTEBOOK_PATH.name}",
        "# This file is a mechanical notebook-section export for refactoring.",
        "# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.",
        "",
    ]
    for cell_index in cells:
        lines.extend(_render_cell(cell_index, notebook["cells"][cell_index]))
    return "\n".join(lines) + "\n"


def _render_readme(unassigned_cells: list[int]) -> str:
    readme = f"""# Notebook Section Exports

These files are a mechanical first pass at converting `../{NOTEBOOK_PATH.name}` into source-controlled Python section files. They are intended to make the existing notebook easier to review and refactor without changing the scientific implementation.

Important constraints:

- The notebook remains the reference implementation.
- Equations, numerical algorithms, parameter values, and figure logic should remain unchanged.
- Code was copied from notebook cells into logical section files; notebook magics such as `%matplotlib inline` are commented out so the exports are valid Python text.
- Some cells appear in more than one section when they currently mix responsibilities, for example camera helpers embedded in imaging/demo cells. Later refactors should deduplicate these mechanically after regression checks.
- These section files are an intermediate refactoring aid, not yet a package API. Execute the original notebook for authoritative results until dependency order has been normalised.

Section files:

"""
    for filename, title, cells in SECTIONS:
        readme += f"- `{filename}` — {title}; notebook cells {list(cells)}.\n"
    readme += f"\nUnassigned notebook cells: {unassigned_cells or 'none'}.\n"
    return readme


def main() -> None:
    notebook = json.loads(NOTEBOOK_PATH.read_text())
    OUTPUT_DIR.mkdir(exist_ok=True)
    assigned_cells: list[int] = []

    for filename, title, cells_iterable in SECTIONS:
        cells = list(cells_iterable)
        assigned_cells.extend(cells)
        (OUTPUT_DIR / filename).write_text(_render_section(filename, title, cells, notebook))

    unassigned_cells = [
        cell_index
        for cell_index in range(len(notebook["cells"]))
        if cell_index not in set(assigned_cells)
    ]
    (OUTPUT_DIR / "README.md").write_text(_render_readme(unassigned_cells))


if __name__ == "__main__":
    main()
