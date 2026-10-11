"""Split a genome by chromosome and compare it against other data.

This ports the project's original exploratory notebook analysis into a tested,
importable module:

* :func:`summarize` — how many markers each chromosome carries, and how many of
  them a trait database (SNPedia / GWAS-style table) actually annotates.
* :func:`compare_genomes` — line up two genomes (e.g. a person and their parent)
  and split every shared SNP into *concordant* (identical call) and *discordant*
  (different call) buckets, per chromosome.

Everything here is pure and dependency-light (no pandas, no OpenCV), so it runs
in the test suite and on a headless box.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .genetics import normalize_genotype
from .genome import iter_snp_records
from .iids import IidResolver
from .models import SnpediaEntry

# Karyotype order: autosomes 1-22, then the sex chromosomes, the pseudoautosomal
# region (AncestryDNA's "XY"/code 25) and mitochondria.
_CHROMOSOME_ORDER = {
    name: i
    for i, name in enumerate(
        [str(n) for n in range(1, 23)] + ["X", "Y", "XY", "MT"]
    )
}


def _chromosome_sort_key(chromosome: str) -> tuple[int, str]:
    """Sort autosomes numerically, then X/Y/MT, then anything unexpected."""
    return (_CHROMOSOME_ORDER.get(chromosome, len(_CHROMOSOME_ORDER)), chromosome)


# --------------------------------------------------------------------------- #
# Per-chromosome summary
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, slots=True)
class ChromosomeStat:
    """Marker tally for one chromosome."""

    chromosome: str
    markers: int      # every row the chip reported on this chromosome
    diploid: int      # usable diploid SNP calls
    matched: int      # of those, ones the trait table annotates


def summarize(
    path: str | Path,
    table: dict[tuple[str, str], SnpediaEntry] | None = None,
    *,
    resolver: IidResolver | None = None,
) -> list[ChromosomeStat]:
    """Split ``path`` by chromosome and tally markers (optionally trait matches).

    When ``table`` is given, every diploid call is looked up in it — resolving
    23andMe internal ids to RSIDs first when ``resolver`` is supplied — so the
    ``matched`` column mirrors the genome × SNPedia join the mirror itself does.
    """
    markers: dict[str, int] = defaultdict(int)
    diploid: dict[str, int] = defaultdict(int)
    matched: dict[str, int] = defaultdict(int)

    for rsid, chromosome, _position, genotype in iter_snp_records(path):
        markers[chromosome] += 1
        key = normalize_genotype(genotype)
        if key is None:
            continue
        diploid[chromosome] += 1
        if table is not None:
            marker = resolver.resolve(rsid) if resolver is not None else rsid
            if (marker, key) in table or (rsid, key) in table:
                matched[chromosome] += 1

    return [
        ChromosomeStat(chrom, markers[chrom], diploid[chrom], matched[chrom])
        for chrom in sorted(markers, key=_chromosome_sort_key)
    ]


def write_summary_csv(stats: list[ChromosomeStat], path: str | Path) -> Path:
    """Write a per-chromosome summary to ``path`` and return it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chromosome", "markers", "diploid", "matched"])
        for s in stats:
            writer.writerow([s.chromosome, s.markers, s.diploid, s.matched])
    return path


# --------------------------------------------------------------------------- #
# Two-genome comparison
# --------------------------------------------------------------------------- #
@dataclass(frozen=True, slots=True)
class DiscordantSnp:
    """A SNP where two genomes report different diploid genotypes."""

    rsid: str
    chromosome: str
    position: str
    genotype_a: str
    genotype_b: str


@dataclass(slots=True)
class ComparisonResult:
    """Outcome of lining up two genomes at their shared SNPs."""

    shared: int = 0                       # SNPs both genomes call as diploid
    concordant: int = 0                   # identical genotype
    discordant_snps: list[DiscordantSnp] = field(default_factory=list)
    #: chromosome -> (shared, concordant, discordant)
    by_chromosome: dict[str, tuple[int, int, int]] = field(default_factory=dict)

    @property
    def discordant(self) -> int:
        return len(self.discordant_snps)

    @property
    def concordance_rate(self) -> float:
        return self.concordant / self.shared if self.shared else 0.0


def _index(path: str | Path) -> dict[str, tuple[str, str, str]]:
    """rsid -> (chromosome, position, normalized_genotype) for diploid calls."""
    index: dict[str, tuple[str, str, str]] = {}
    for rsid, chromosome, position, genotype in iter_snp_records(path):
        key = normalize_genotype(genotype)
        if key is None:
            continue
        index.setdefault(rsid, (chromosome, position, key))
    return index


def compare_genomes(
    path_a: str | Path, path_b: str | Path
) -> ComparisonResult:
    """Compare two genomes SNP-by-SNP, bucketing shared calls per chromosome.

    Only SNPs both files call as clean diploids are considered.  A *discordant*
    SNP is one present in both with different genotypes — the raw material for
    parent/child or duplicate-sample checks (as in the notebook's ``lineage``
    pass, minus the IBD-segment reconstruction).
    """
    index_a = _index(path_a)
    index_b = _index(path_b)

    result = ComparisonResult()
    per_chrom_shared: dict[str, int] = defaultdict(int)
    per_chrom_conc: dict[str, int] = defaultdict(int)
    per_chrom_disc: dict[str, int] = defaultdict(int)

    for rsid, (chromosome, position, geno_a) in index_a.items():
        other = index_b.get(rsid)
        if other is None:
            continue
        _chrom_b, _pos_b, geno_b = other
        result.shared += 1
        per_chrom_shared[chromosome] += 1
        if geno_a == geno_b:
            result.concordant += 1
            per_chrom_conc[chromosome] += 1
        else:
            per_chrom_disc[chromosome] += 1
            result.discordant_snps.append(
                DiscordantSnp(rsid, chromosome, position, geno_a, geno_b)
            )

    result.discordant_snps.sort(
        key=lambda d: (_chromosome_sort_key(d.chromosome), d.rsid)
    )
    result.by_chromosome = {
        chrom: (
            per_chrom_shared[chrom],
            per_chrom_conc[chrom],
            per_chrom_disc[chrom],
        )
        for chrom in sorted(per_chrom_shared, key=_chromosome_sort_key)
    }
    return result


def write_discordant_csv(result: ComparisonResult, path: str | Path) -> Path:
    """Write discordant SNPs to ``path`` and return it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rsid", "chromosome", "position", "genotype_a", "genotype_b"])
        for d in result.discordant_snps:
            writer.writerow(
                [d.rsid, d.chromosome, d.position, d.genotype_a, d.genotype_b]
            )
    return path
