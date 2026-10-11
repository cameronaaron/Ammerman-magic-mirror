"""Fetch gene names for rsids by scraping SNPedia SNP-page wikitext.

SNPedia does not expose ``Gene`` as a queryable Semantic-MediaWiki property, but
each SNP page's raw wikitext carries it in the ``{{Rsnum|…|Gene=APOE|…}}``
template.  This fetches that content for a set of rsids (batched, ~50 titles per
request) and caches ``rsid → gene`` on disk so repeat runs are instant.

Only the *matched* rsids are ever fetched — a few hundred, not all of SNPedia —
so gene enrichment stays cheap.  Any network failure degrades gracefully to
"no genes" rather than raising.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable, Iterable
from pathlib import Path

import httpx

from .live import API_ENDPOINT, USER_AGENT

_GENE_RE = re.compile(r"\bGene\s*=\s*([A-Za-z0-9][A-Za-z0-9._@-]*)")
_BATCH = 50
_CACHE_NAME = "snpedia_genes.json"

Logger = Callable[[str], None]


def _title(rsid: str) -> str:
    return rsid[:1].upper() + rsid[1:]


def _read_cache(path: Path) -> dict[str, str]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _write_cache(path: Path, data: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")


def _extract_gene(wikitext: str) -> str:
    match = _GENE_RE.search(wikitext)
    return match.group(1).strip() if match else ""


def fetch_genes(
    rsids: Iterable[str],
    *,
    cache_dir: Path,
    timeout: float = 60.0,
    log: Logger = print,
) -> dict[str, str]:
    """Return ``rsid → gene`` for ``rsids`` (only real rs-ids), disk-cached.

    Unknown/failed lookups are cached as ``""`` so we don't re-fetch them.  A
    network error returns whatever is already cached rather than raising.
    """
    cache_path = Path(cache_dir) / _CACHE_NAME
    cache = _read_cache(cache_path)

    wanted = [r for r in {r.lower() for r in rsids} if r.startswith("rs")]
    missing = [r for r in wanted if r not in cache]
    if not missing:
        return {r: cache[r] for r in wanted if cache.get(r)}

    log(f"Fetching gene names for {len(missing)} rsids from SNPedia…")
    headers = {"User-Agent": USER_AGENT}
    try:
        with httpx.Client(timeout=timeout, headers=headers, follow_redirects=True) as c:
            for start in range(0, len(missing), _BATCH):
                batch = missing[start : start + _BATCH]
                params = {
                    "action": "query",
                    "format": "json",
                    "formatversion": "2",
                    "prop": "revisions",
                    "rvprop": "content",
                    "rvslots": "main",
                    "titles": "|".join(_title(r) for r in batch),
                }
                data = c.get(API_ENDPOINT, params=params).raise_for_status().json()
                found = _parse_pages(data)
                for rsid in batch:  # cache "" for the ones with no gene / no page
                    cache[rsid] = found.get(rsid, "")
    except (httpx.HTTPError, json.JSONDecodeError) as exc:
        log(f"  ! gene fetch failed ({exc}); using whatever is cached")

    _write_cache(cache_path, cache)
    return {r: cache[r] for r in wanted if cache.get(r)}


def _parse_pages(data: dict) -> dict[str, str]:
    genes: dict[str, str] = {}
    for page in data.get("query", {}).get("pages", []):
        title = page.get("title", "")
        rsid = title.lower()
        if not rsid.startswith("rs"):
            continue
        revisions = page.get("revisions", [])
        if not revisions:
            continue
        content = revisions[0].get("slots", {}).get("main", {}).get("content", "")
        gene = _extract_gene(content)
        if gene:
            genes[rsid] = gene
    return genes
