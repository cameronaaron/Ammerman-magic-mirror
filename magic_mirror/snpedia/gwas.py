"""GWAS Catalog association enrichment — real odds ratios and risk-allele
population frequencies for a set of rsids.

Promethease's second report layers genome-wide-association data (which allele
raises risk for a trait, and how common that allele is) on top of the per-SNP
summaries.  This fetches that from the EBI **GWAS Catalog** REST API for a
bounded set of rsids (your notable SNPs), caches it, and lets the caller work
out how many risk alleles you carry.

It is deliberately scoped to the rsids handed in (not a full-genome scan, which
would need the multi-hundred-MB catalog download) and degrades gracefully to
"no data" on any network error.
"""

from __future__ import annotations

import json
import re
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

import httpx

from .live import USER_AGENT

GWAS_API = (
    "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/"
    "{rsid}/associations?projection=associationBySnp"
)
_ALLELE_RE = re.compile(r"-([ACGT])$")
_CACHE_NAME = "gwas_associations.json"

Logger = Callable[[str], None]


@dataclass(frozen=True, slots=True)
class GwasAssoc:
    """One GWAS Catalog association for a SNP."""

    rsid: str
    trait: str
    risk_allele: str          # single base
    risk_freq: float | None   # population frequency of the risk allele (RAF)
    pvalue: float | None
    effect: str               # OR-per-copy or beta, as reported

    def dosage(self, genotype: str) -> int:
        """How many copies of the risk allele a normalized genotype carries."""
        return genotype.upper().count(self.risk_allele)


def _read_cache(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _to_float(value: object) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _parse_response(rsid: str, data: dict) -> list[list]:
    """Flatten one SNP's associations into cache rows (best p-value per trait)."""
    best: dict[tuple[str, str], list] = {}
    for assoc in data.get("_embedded", {}).get("associations", []):
        pvalue = _to_float(assoc.get("pvalue"))
        raf = _to_float(assoc.get("riskFrequency"))
        effect = assoc.get("orPerCopyNum") or assoc.get("betaNum") or ""
        traits = [t.get("trait", "") for t in assoc.get("efoTraits", [])]
        alleles = [
            ra.get("riskAlleleName", "")
            for locus in assoc.get("loci", [])
            for ra in locus.get("strongestRiskAlleles", [])
        ]
        for name in alleles:
            m = _ALLELE_RE.search(name or "")
            if not m:
                continue
            allele = m.group(1)
            for trait in traits:
                if not trait:
                    continue
                key = (trait, allele)
                prev = best.get(key)
                stronger = pvalue is not None and (
                    prev is None or prev[3] is None or pvalue < prev[3]
                )
                if prev is None or stronger:
                    best[key] = [trait, allele, raf, pvalue, str(effect)]
    return list(best.values())


def fetch_gwas(
    rsids: Iterable[str],
    *,
    cache_dir: Path,
    timeout: float = 60.0,
    log: Logger = print,
) -> dict[str, list[GwasAssoc]]:
    """Return ``rsid → [GwasAssoc]`` for real rs-ids, disk-cached per rsid."""
    cache_path = Path(cache_dir) / _CACHE_NAME
    cache = _read_cache(cache_path)

    wanted = [r for r in {r.lower() for r in rsids} if r.startswith("rs")]
    missing = [r for r in wanted if r not in cache]
    if missing:
        log(f"Fetching GWAS Catalog associations for {len(missing)} rsids…")
        headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
        try:
            with httpx.Client(
                timeout=timeout, headers=headers, follow_redirects=True
            ) as client:
                for rsid in missing:
                    try:
                        resp = client.get(GWAS_API.format(rsid=rsid))
                        if resp.status_code == 404:
                            cache[rsid] = []
                            continue
                        resp.raise_for_status()
                        cache[rsid] = _parse_response(rsid, resp.json())
                    except (httpx.HTTPError, json.JSONDecodeError):
                        cache[rsid] = []  # skip this one, keep going
                    time.sleep(0.1)  # be gentle on the EBI API
        finally:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache), encoding="utf-8")

    return {
        rsid: [GwasAssoc(rsid, *row) for row in cache.get(rsid, [])]
        for rsid in wanted
        if cache.get(rsid)
    }
