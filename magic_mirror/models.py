"""Domain value objects shared across the pipeline.

Plain frozen dataclasses (not pydantic) — these sit on the hot path and carry no
validation cost; pydantic is reserved for user-facing config in :mod:`config`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class GenomeData:
    """A parsed genome: rsid -> normalized genotype, plus provenance counts."""

    calls: dict[str, str]          # rsid (lower) -> normalized genotype
    total_markers: int             # every marker the chip reported
    diploid_markers: int           # usable diploid SNP calls
    source: Path
    internal_ids: int = 0          # markers reported under 23andMe "i" ids
    resolved_ids: int = 0          # of those, translated to a standard RSID

    def __len__(self) -> int:
        return self.diploid_markers


@dataclass(frozen=True, slots=True)
class SnpediaEntry:
    """One genotype SNPedia documents: magnitude, repute and a summary."""

    rsid: str          # lower-case
    genotype: str      # normalized (sorted alleles)
    magnitude: float
    summary: str
    repute: str = ""   # SNPedia "Repute": "Good", "Bad" or "" (unset/neutral)

    @property
    def key(self) -> tuple[str, str]:
        return (self.rsid, self.genotype)


@dataclass(frozen=True, slots=True)
class Trait:
    """A genotype the user carries that SNPedia has something to say about."""

    rsid: str
    genotype: str
    magnitude: float
    summary: str
    category: str
    repute: str = "Neutral"  # "Good" | "Bad" | "Neutral"

    @property
    def display_genotype(self) -> str:
        """The genotype as a reader expects it, e.g. ``"G;G"``."""
        return ";".join(self.genotype)

    @property
    def snpedia_url(self) -> str:
        return f"https://www.snpedia.com/index.php/{self.rsid}"

    @property
    def genotype_url(self) -> str:
        """Link to the specific genotype page, e.g. ``Rs53576(G;G)``."""
        head = self.rsid[:1].upper() + self.rsid[1:]
        return f"https://www.snpedia.com/index.php/{head}({self.display_genotype})"


@dataclass(frozen=True, slots=True)
class FaceBox:
    """A detected face rectangle in frame pixel coordinates."""

    x: int
    y: int
    w: int
    h: int

    @property
    def center(self) -> tuple[int, int]:
        return self.x + self.w // 2, self.y + self.h // 2
