from __future__ import annotations

from pathlib import Path
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
DELIVERABLES_DIR = REPO_ROOT / "deliverables"
BUNDLE_PATH = DELIVERABLES_DIR / "non_destructive_image_code_bundle.zip"

INCLUDE_PATHS = [
    ".gitignore",
    "README.md",
    "requirements.txt",
    "pytest.ini",
    "1 calculations revised 2  multishot  6  extended.ipynb",
    "deliverables/CODE_BUNDLE_README.md",
    "docs",
    "notebook_sections",
    "src",
    "scripts",
    "tests",
    "regression",
]

EXCLUDED_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".idea",
    ".vscode",
    "node_modules",
}

EXCLUDED_FILENAMES = {
    ".DS_Store",
    "Thumbs.db",
}

EXCLUDED_SUFFIXES = {
    ".pyc",
    ".pyo",
    ".tmp",
    ".swp",
    ".swo",
}

ZIP_TIMESTAMP = (2026, 1, 1, 0, 0, 0)


def should_exclude(path: Path) -> bool:
    parts = set(path.parts)
    if parts & EXCLUDED_DIRS:
        return True
    if path.name in EXCLUDED_FILENAMES:
        return True
    if path.suffix in EXCLUDED_SUFFIXES:
        return True
    if path.parts[:1] == ("deliverables",) and path.suffix == ".zip":
        return True
    if path.parts[:2] == ("regression", "baseline") and path.suffix == ".npz":
        return True
    return False


def iter_bundle_files() -> list[Path]:
    files: list[Path] = []
    for include_path in INCLUDE_PATHS:
        root = REPO_ROOT / include_path
        if not root.exists():
            raise FileNotFoundError(f"Required bundle path is missing: {include_path}")
        if root.is_file():
            relative = root.relative_to(REPO_ROOT)
            if not should_exclude(relative):
                files.append(root)
            continue
        for path in root.rglob("*"):
            if path.is_file():
                relative = path.relative_to(REPO_ROOT)
                if not should_exclude(relative):
                    files.append(path)
    return sorted(files, key=lambda path: path.relative_to(REPO_ROOT).as_posix())


def read_reproducible_bytes(path: Path) -> bytes:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return data
    return text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8")


def write_bundle() -> Path:
    DELIVERABLES_DIR.mkdir(exist_ok=True)
    if BUNDLE_PATH.exists():
        BUNDLE_PATH.unlink()

    with zipfile.ZipFile(BUNDLE_PATH, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in iter_bundle_files():
            relative = path.relative_to(REPO_ROOT).as_posix()
            info = zipfile.ZipInfo(relative, ZIP_TIMESTAMP)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, read_reproducible_bytes(path))

    return BUNDLE_PATH


def main() -> None:
    bundle = write_bundle()
    print(bundle.relative_to(REPO_ROOT).as_posix())


if __name__ == "__main__":
    main()

