"""Validate the mechanical notebook section exports.

The validation is intentionally lightweight and non-scientific: for unmigrated
sections it confirms that `notebook_sections/` can be regenerated exactly from
the notebook JSON, and for migrated sections it still checks Python syntax. It
does not execute the simulation or compare numerical physics outputs.
"""

from __future__ import annotations

import importlib.util
import json
import py_compile
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
EXPORT_SCRIPT = REPO_ROOT / "scripts" / "export_notebook_sections.py"

# Sections intentionally migrated away from exact mechanical export text. These
# remain behaviour-preserving targets but no longer compare byte-for-byte to the
# original notebook cell source.
MIGRATED_SECTIONS = {"02_atomic_model.py", "03_light_atom_interaction.py", "04_pci.py", "06_faraday.py"}


def _load_export_module():
    spec = importlib.util.spec_from_file_location("export_notebook_sections", EXPORT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {EXPORT_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    exporter = _load_export_module()
    notebook = json.loads((REPO_ROOT / exporter.NOTEBOOK_PATH).read_text())
    errors: list[str] = []
    assigned_cells: list[int] = []

    for filename, title, cells_iterable in exporter.SECTIONS:
        cells = list(cells_iterable)
        assigned_cells.extend(cells)
        expected = exporter._render_section(filename, title, cells, notebook)
        section_path = REPO_ROOT / exporter.OUTPUT_DIR / filename
        if not section_path.exists():
            errors.append(f"Missing exported section: {section_path}")
            continue
        actual = section_path.read_text()
        if filename not in MIGRATED_SECTIONS and actual != expected:
            errors.append(f"Out-of-date exported section: {section_path}")
        try:
            py_compile.compile(str(section_path), doraise=True)
        except py_compile.PyCompileError as exc:
            errors.append(f"Python syntax error in {section_path}: {exc.msg}")

    unassigned_cells = [
        cell_index
        for cell_index in range(len(notebook["cells"]))
        if cell_index not in set(assigned_cells)
    ]
    expected_readme = exporter._render_readme(unassigned_cells)
    readme_path = REPO_ROOT / exporter.OUTPUT_DIR / "README.md"
    if readme_path.read_text() != expected_readme:
        errors.append(f"Out-of-date exported section README: {readme_path}")

    if errors:
        print("Notebook section validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    migrated = ", ".join(sorted(MIGRATED_SECTIONS)) or "none"
    print(
        "Notebook section exports are in sync for unmigrated sections; "
        f"migrated sections ({migrated}) are syntactically valid."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
