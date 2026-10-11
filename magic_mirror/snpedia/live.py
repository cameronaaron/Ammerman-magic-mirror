"""Live SNPedia client — pull the notable-genotype table straight from the wiki.

SNPedia runs Semantic MediaWiki, so we can ask for every genotype page above a
magnitude threshold along with its ``Magnitude`` and ``Summary`` properties.
That scan is authoritative but slow (~10 s per 500-row page), so results are
cached to disk; the first run populates the cache and later runs are instant.

HTTP is via ``httpx`` with a bounded retry loop for the wiki's occasional 5xx.
"""

from __future__ import annotations

import json
import time
from collections.abc import Callable
from pathlib import Path

import httpx

from ..genetics import coerce_magnitude, parse_genotype_title
from ..models import SnpediaEntry
from .parse import Table

API_ENDPOINT = "https://bots.snpedia.com/api.php"
USER_AGENT = (
    "GeneticReflexions/2.0 (magic-mirror; https://github.com/; local research use)"
)
_PAGE_LIMIT = 500
_MAX_RETRIES = 3
_RETRY_BACKOFF = 2.0
_INTER_PAGE_PAUSE = 0.5  # be polite to a volunteer-run wiki
_STALL_LIMIT = 2         # consecutive no-progress pages before giving up

Logger = Callable[[str], None]


class SnpediaUnavailable(RuntimeError):
    """Raised when the live SNPedia API cannot be reached."""


# --------------------------------------------------------------------------- #
# Cache
# --------------------------------------------------------------------------- #
def _cache_path(cache_dir: Path, min_magnitude: float) -> Path:
    return cache_dir / f"snpedia_mag_{min_magnitude:g}.json"


def _read_cache(path: Path, ttl_seconds: float) -> list | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if time.time() - payload.get("fetched_at", 0) > ttl_seconds:
        return None
    return payload.get("entries")


def _write_cache(path: Path, min_magnitude: float, entries: list[SnpediaEntry]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "fetched_at": time.time(),
        "min_magnitude": min_magnitude,
        "source": API_ENDPOINT,
        "entries": [
            [e.rsid, e.genotype, e.magnitude, e.summary, e.repute] for e in entries
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _rows_to_table(rows: list) -> Table:
    table: Table = {}
    for row in rows:
        rsid, genotype, magnitude, summary = row[:4]
        repute = row[4] if len(row) > 4 else ""  # tolerate pre-Repute caches
        entry = SnpediaEntry(rsid, genotype, float(magnitude), summary, repute)
        table[entry.key] = entry
    return table


# --------------------------------------------------------------------------- #
# HTTP
# --------------------------------------------------------------------------- #
def _ask(client: httpx.Client, query: str, offset: int) -> dict:
    params = {"action": "ask", "format": "json", "query": f"{query}|offset={offset}"}
    last_error: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            response = client.get(API_ENDPOINT, params=params)
            response.raise_for_status()
            return response.json()
        except (httpx.HTTPError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < _MAX_RETRIES - 1:
                time.sleep(_RETRY_BACKOFF * (attempt + 1))
    raise SnpediaUnavailable(f"SNPedia API request failed: {last_error}")


def _walk(query: str, timeout: float, log: Logger) -> list[SnpediaEntry]:
    """Paginate a Semantic MediaWiki ``ask`` query into ``SnpediaEntry`` rows.

    Shared by the notable-genotype fetch and the full-snapshot walk; a tail
    failure keeps whatever was already gathered rather than discarding it.
    """
    entries: dict[tuple[str, str], SnpediaEntry] = {}
    offset = 0
    stalled = 0
    headers = {"User-Agent": USER_AGENT}
    with httpx.Client(timeout=timeout, headers=headers, follow_redirects=True) as client:
        while True:
            try:
                data = _ask(client, query, offset)
            except SnpediaUnavailable as exc:
                # A tail failure shouldn't discard a whole successful scan.
                if entries:
                    log(f"  ! stopped early ({exc}); keeping {len(entries)} genotypes")
                    break
                raise

            results = data.get("query", {}).get("results", {})
            if not results:
                break

            before = len(entries)
            for title, body in results.items():
                parsed = parse_genotype_title(title)
                if parsed is None:
                    continue
                rsid, genotype = parsed
                printouts = body.get("printouts", {})
                magnitude = coerce_magnitude(printouts.get("Magnitude"))
                summary = _first_str(printouts.get("Summary"))
                repute = _first_str(printouts.get("Repute"))
                entries[(rsid, genotype)] = SnpediaEntry(
                    rsid, genotype, magnitude, summary, repute
                )
            log(f"  …{len(entries)} genotypes so far")

            # Guard against a continuation offset that plateaus (would loop forever).
            stalled = stalled + 1 if len(entries) == before else 0
            if stalled >= _STALL_LIMIT:
                break
            cont = data.get("query-continue-offset")
            if not cont or int(cont) <= offset:
                break
            offset = int(cont)
            time.sleep(_INTER_PAGE_PAUSE)
    return list(entries.values())


def _download(min_magnitude: float, timeout: float, log: Logger) -> list[SnpediaEntry]:
    # `>` is strict, so nudge the floor down slightly to keep entries exactly at
    # the threshold; the display filter re-applies the precise cutoff later.
    floor = max(0.0, min_magnitude - 0.001)
    query = (
        f"[[Category:Is a genotype]][[Magnitude::>{floor:g}]]"
        f"|?Magnitude|?Repute|?Summary|limit={_PAGE_LIMIT}"
    )
    log(f"Fetching notable genotypes from SNPedia (magnitude ≥ {min_magnitude:g})…")
    entries = _walk(query, timeout, log)
    if not entries:
        raise SnpediaUnavailable("SNPedia returned no genotypes for that threshold.")
    return entries


# Semantic MediaWiki caps how far a single ``ask`` query paginates (~5000 rows),
# so a naive "all genotypes" walk silently truncates alphabetically.  Walking in
# magnitude bands keeps every band under the cap; the union is comprehensive.
_SNAPSHOT_BANDS = [-1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 100]


def fetch_snapshot(*, timeout: float = 60.0, log: Logger = print) -> list[SnpediaEntry]:
    """Walk SNPedia genotypes in magnitude bands for a comprehensive snapshot.

    Returns every genotype that carries a human-readable summary (the only ones
    the app can use), with Magnitude and Repute — the data behind ``result.csv``.
    """
    log("Walking SNPedia genotypes in magnitude bands for a fresh snapshot…")
    merged: dict[tuple[str, str], SnpediaEntry] = {}
    for lo, hi in zip(_SNAPSHOT_BANDS, _SNAPSHOT_BANDS[1:], strict=False):
        # Tiny overlap + dict-dedup so a genotype exactly on a band edge is never
        # lost between (lo, hi) intervals.
        query = (
            f"[[Category:Is a genotype]][[Magnitude::>{lo - 0.001:g}]]"
            f"[[Magnitude::<{hi + 0.001:g}]]"
            f"|?Magnitude|?Repute|?Summary|limit={_PAGE_LIMIT}"
        )
        for entry in _walk(query, timeout, log):
            merged[entry.key] = entry
        log(f"  band ({lo:g}–{hi:g}]: {len(merged)} genotypes so far")

    # Keep only genotypes with a summary — build_traits requires one, so
    # summary-less rows would just bloat the snapshot (as the old one did).
    entries = [e for e in merged.values() if e.summary]
    if not entries:
        raise SnpediaUnavailable("SNPedia returned no genotypes.")
    return entries


def _first_str(value: object) -> str:
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    return str(value or "").strip()


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def fetch_live(
    min_magnitude: float,
    *,
    cache_dir: Path,
    ttl_seconds: float,
    timeout: float = 60.0,
    force_refresh: bool = False,
    log: Logger = print,
) -> Table:
    """Return the live SNPedia trait table, using an on-disk cache when fresh.

    Raises:
        SnpediaUnavailable: if the API is unreachable and no cache exists.
    """
    path = _cache_path(Path(cache_dir), min_magnitude)

    if not force_refresh:
        cached = _read_cache(path, ttl_seconds)
        if cached is not None:
            log(f"Using cached SNPedia data ({len(cached)} genotypes).")
            return _rows_to_table(cached)

    entries = _download(min_magnitude, timeout, log)
    _write_cache(path, min_magnitude, entries)
    log(f"Cached {len(entries)} genotypes → {path}")
    return {e.key: e for e in entries}
