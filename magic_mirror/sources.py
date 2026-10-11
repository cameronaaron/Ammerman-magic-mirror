"""Discover a genome file on inserted removable media.

The installation's intended UX is "insert a USB drive containing your genome
file".  This module locates that file across the usual mount points on macOS
(``/Volumes``) and Linux (``/media/<user>``, ``/run/media/<user>``), verifying it
really looks like a 23andMe or AncestryDNA export before offering it.
"""

from __future__ import annotations

import getpass
from collections.abc import Iterator
from pathlib import Path

# Filenames 23andMe / AncestryDNA / common tools produce, in rough priority order.
_GENOME_GLOBS = (
    "*genome*.txt", "*23andme*.txt", "*ancestry*.txt", "genome*.txt", "*.txt"
)
_MAX_SCAN_BYTES = 4096  # only sniff the header to validate


def _mount_roots() -> Iterator[Path]:
    user = getpass.getuser()
    candidates = [
        Path("/Volumes"),                 # macOS
        Path("/media") / user,            # Linux (Ubuntu/Debian)
        Path("/run/media") / user,        # Linux (Fedora/Arch)
        Path("/media"),                   # Linux (fallback)
    ]
    for root in candidates:
        if root.is_dir():
            yield root


def looks_like_23andme(path: Path) -> bool:
    """Cheap header sniff: is this plausibly a 23andMe or AncestryDNA raw file?

    (Name kept for backwards compatibility; it now accepts either vendor.)
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            head = handle.read(_MAX_SCAN_BYTES).lower()
    except OSError:
        return False
    # A vendor banner, or the canonical column header of either format.
    if "23andme" in head or "ancestrydna" in head:
        return True
    return "rsid" in head and ("genotype" in head or "allele1" in head)


def iter_candidate_genomes() -> Iterator[Path]:
    """Yield genome-like files found on mounted removable media."""
    seen: set[Path] = set()
    for root in _mount_roots():
        for volume in sorted(root.iterdir()) if root.is_dir() else []:
            if not volume.is_dir():
                continue
            for pattern in _GENOME_GLOBS:
                for match in sorted(volume.glob(pattern)):
                    if match in seen or not match.is_file():
                        continue
                    seen.add(match)
                    if looks_like_23andme(match):
                        yield match


def find_usb_genome() -> Path | None:
    """Return the first plausible genome on removable media, or ``None``."""
    return next(iter_candidate_genomes(), None)


def iter_folder_genomes(folder: str | Path) -> Iterator[Path]:
    """Yield genome-like files dropped in a local folder (non-recursive).

    The installation mostly reads from USB, but a local drop folder in the
    project directory is handier for testing and for a fixed kiosk setup.
    """
    folder = Path(folder)
    if not folder.is_dir():
        return
    seen: set[Path] = set()
    for pattern in _GENOME_GLOBS:
        for match in sorted(folder.glob(pattern)):
            if match in seen or not match.is_file():
                continue
            seen.add(match)
            if looks_like_23andme(match):
                yield match


def find_local_genome(folder: str | Path) -> Path | None:
    """Return the first plausible genome in ``folder``, or ``None``."""
    return next(iter_folder_genomes(folder), None)
