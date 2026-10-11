"""Load and write the frozen SNPedia snapshot (``result.csv``)."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

from ..genetics import coerce_magnitude, parse_genotype_title
from ..models import SnpediaEntry

Table = dict[tuple[str, str], SnpediaEntry]


def _title(entry: SnpediaEntry) -> str:
    """Reconstruct a SNPedia genotype-page title, e.g. ``Rs17822931(T;T)``."""
    rsid = entry.rsid[:1].upper() + entry.rsid[1:]
    alleles = ";".join(entry.genotype)
    return f"{rsid}({alleles})"


def write_snapshot_csv(entries: Iterable[SnpediaEntry], path: str | Path) -> Path:
    """Write a SNPedia snapshot to ``result.csv`` form, now with a Repute column.

    Header is ``,Magnitude,Repute,Summary``; :func:`load_csv` reads the Repute
    column back automatically.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted(entries, key=lambda e: (-e.magnitude, e.rsid, e.genotype))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", "Magnitude", "Repute", "Summary"])
        for e in rows:
            writer.writerow([_title(e), f"{e.magnitude:g}", e.repute, e.summary])
    return path


def load_csv(path: str | Path) -> Table:
    """Return a ``(rsid, genotype) -> SnpediaEntry`` table from ``result.csv``.

    Rows whose key isn't a clean ``rsNNN(X;Y)`` / ``iNNN(X;Y)`` diploid pattern
    (indels, malformed rows) are skipped rather than guessed at.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"SNPedia table not found: {path}")

    table: Table = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle)
        # Columns are found by name: snapshots exist as ",Magnitude,Summary"
        # and as ",Magnitude,Repute,Summary" (what write_snapshot_csv writes).
        header = next(reader, None) or []
        lower = [h.strip().lower() for h in header]
        repute_idx = lower.index("repute") if "repute" in lower else None
        summary_idx = lower.index("summary") if "summary" in lower else 2
        for row in reader:
            if not row or not row[0]:
                continue
            parsed = parse_genotype_title(row[0])
            if parsed is None:
                continue
            rsid, genotype = parsed
            magnitude = coerce_magnitude(row[1] if len(row) > 1 else "")
            summary = (row[summary_idx] if len(row) > summary_idx else "").strip()
            repute = ""
            if repute_idx is not None and len(row) > repute_idx:
                repute = row[repute_idx].strip()
            table[(rsid, genotype)] = SnpediaEntry(
                rsid, genotype, magnitude, summary, repute
            )

    if not table:
        raise ValueError(f"{path} yielded no usable SNPedia rows.")
    return table
