"""Capture the executable source and runtime state of reconstruction studies."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import platform
from typing import Iterable

from .io import file_sha256, git_value


_SHARED_SOURCE_FILES = (
    Path("pyproject.toml"),
    Path("src/non_destructive_image/camera.py"),
    Path("src/non_destructive_image/fourier.py"),
)


def _package_version(distribution: str) -> str:
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "not-installed"


def reconstruction_source_hashes(
    repository_root: Path,
    *,
    entry_points: Iterable[Path] = (),
) -> dict[str, str]:
    """Hash every in-repository Python source used by a reconstruction study."""

    repository_root = repository_root.resolve()
    package_root = repository_root / "src/non_destructive_image"
    paths = {
        *(repository_root / relative for relative in _SHARED_SOURCE_FILES),
        *(repository_root / relative for relative in entry_points),
        *package_root.rglob("*.py"),
    }
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "provenance source file is missing: " + ", ".join(map(str, missing))
        )
    return {
        path.relative_to(repository_root).as_posix(): file_sha256(path)
        for path in sorted(paths, key=lambda item: item.as_posix())
    }


def capture_reconstruction_provenance(
    repository_root: Path,
    *,
    entry_points: Iterable[Path] = (),
) -> dict[str, object]:
    """Record Git disclosure, source hashes and numerical-library versions."""

    repository_root = repository_root.resolve()
    git_status = git_value(repository_root, "status", "--porcelain")
    return {
        "git_branch_before_write": git_value(
            repository_root, "branch", "--show-current"
        ),
        "git_commit_before_write": git_value(repository_root, "rev-parse", "HEAD"),
        "git_status_porcelain_before_write": git_status,
        "git_dirty_before_write": bool(git_status),
        "runtime_versions": {
            "python": platform.python_version(),
            "numpy": _package_version("numpy"),
            "scipy": _package_version("scipy"),
            "matplotlib": _package_version("matplotlib"),
        },
        "source_files_sha256": reconstruction_source_hashes(
            repository_root,
            entry_points=entry_points,
        ),
    }
