"""Genetic Reflexions — a magic mirror that displays your genome, not your face.

Layered so each piece is useful (and testable) on its own:

- :mod:`magic_mirror.genetics`  — pure genotype helpers.
- :mod:`magic_mirror.genome`    — parse a 23andMe raw-data file.
- :mod:`magic_mirror.snpedia`   — trait tables (live API or CSV snapshot).
- :mod:`magic_mirror.traits`    — join genome × SNPedia into ranked traits.
- :mod:`magic_mirror.sources`   — discover a genome on removable media.
- :mod:`magic_mirror.theme` / :mod:`~magic_mirror.render` — Pillow compositor.
- :mod:`magic_mirror.detect` / :mod:`~magic_mirror.camera` — optional OpenCV.
- :mod:`magic_mirror.app`       — orchestration (live mirror + still poster).
- :mod:`magic_mirror.cli`       — Typer + Rich command line.

Only the live-mirror path needs OpenCV; everything else runs on Pillow + numpy.
"""

from __future__ import annotations

from .config import Settings, Source
from .models import FaceBox, GenomeData, SnpediaEntry, Trait

__all__ = [
    "FaceBox",
    "GenomeData",
    "Settings",
    "SnpediaEntry",
    "Source",
    "Trait",
    "__version__",
]

__version__ = "2.0.0"
