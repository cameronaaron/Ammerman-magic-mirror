"""SNPedia trait tables from either the live API or the frozen CSV snapshot."""

from __future__ import annotations

from ..config import Settings, Source
from .live import SnpediaUnavailable, fetch_live, fetch_snapshot
from .parse import Table, load_csv, write_snapshot_csv

__all__ = [
    "SnpediaUnavailable",
    "Table",
    "fetch_live",
    "fetch_snapshot",
    "load_csv",
    "resolve_table",
    "write_snapshot_csv",
]


def resolve_table(settings: Settings, *, log=print) -> Table:
    """Resolve the SNPedia trait table per ``settings.source``.

    "live"/"auto" query SNPedia (disk-cached) and fall back to the CSV snapshot
    when the network is unavailable; "csv" uses the snapshot directly.
    """
    if settings.source is Source.CSV:
        return load_csv(settings.snpedia_path)

    # Fetch at a fixed floor (independent of the display threshold) so tweaking
    # --min-magnitude reuses the cache; only a threshold *below* the floor forces
    # a wider scan.
    floor = min(settings.min_magnitude, settings.snpedia_fetch_floor)
    try:
        return fetch_live(
            floor,
            cache_dir=settings.cache_dir,
            ttl_seconds=settings.cache_ttl_seconds,
            timeout=settings.http_timeout,
            force_refresh=settings.force_refresh,
            log=log,
        )
    except SnpediaUnavailable:
        # "live" only falls back when a snapshot exists; otherwise surface it.
        if settings.source is Source.LIVE and not settings.snpedia_path.exists():
            raise
        log("! SNPedia live fetch failed; using local snapshot.")
        return load_csv(settings.snpedia_path)
