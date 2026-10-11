"""Parse a consumer DNA raw-data export into genotype calls.

Two tab-separated formats are supported, told apart by column count per row:

* **23andMe** — ``rsid  chromosome  position  genotype`` (4 columns; the
  genotype is a two-character string like ``"AG"``), chromosomes as
  ``1..22, X, Y, MT``.
* **AncestryDNA** — ``rsid  chromosome  position  allele1  allele2`` (5 columns;
  the two alleles are separate), chromosomes as numeric codes
  (``23=X, 24=Y, 25=XY, 26=MT``) and ``0`` alleles for no-calls.

Both carry a ``# ...`` comment header. We only need rsid -> genotype, so the
parser streams the file into a compact dict rather than a full DataFrame — the
mirror does dictionary lookups, not table scans, so this is both faster and
dependency-free for 600k+ rows.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from .genetics import normalize_genotype
from .iids import IidResolver, is_internal_id
from .models import GenomeData

try:  # optional native (Rust) parser — see rust/crates/genome-py
    import genome_core as _native
except ImportError:  # pragma: no cover - pure-Python fallback path
    _native = None

# AncestryDNA encodes the sex chromosomes and mitochondria numerically; map them
# to the same labels 23andMe uses so downstream code sees one vocabulary.
_ANCESTRY_CHROM = {"23": "X", "24": "Y", "25": "XY", "26": "MT"}


def _iter_data_lines(path: Path) -> Iterator[list[str]]:
    """Yield normalised ``[rsid, chromosome, position, genotype]`` rows.

    Handles both the 23andMe (4-column) and AncestryDNA (5-column) layouts,
    combining Ancestry's split alleles and remapping its numeric chromosome
    codes so the rest of the pipeline is format-agnostic.
    """
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            parts = [p.strip() for p in raw.rstrip("\n").split("\t")]
            # Header row of either format ("rsid ..."); 23andMe even prefixes it
            # with a stray space, already removed by the strip above.
            if not parts or parts[0].lower() == "rsid":
                continue
            if len(parts) == 4:                     # 23andMe
                rsid, chrom, pos, genotype = parts
            elif len(parts) == 5:                   # AncestryDNA
                rsid, chrom, pos, allele1, allele2 = parts
                genotype = allele1 + allele2
                chrom = _ANCESTRY_CHROM.get(chrom, chrom)
            else:
                continue
            yield [rsid, chrom, pos, genotype]


def iter_snp_records(path: str | Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield ``(rsid, chromosome, position, genotype)`` rows from a genome file.

    Accepts either 23andMe or AncestryDNA exports.  Unlike :func:`load_genome`,
    this preserves the chromosome and position columns — needed by the
    per-chromosome tooling — and does not normalise or filter genotypes.  The
    header and comment lines are skipped.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Genome file not found: {path}")
    for rsid, chromosome, position, genotype in _iter_data_lines(path):
        yield rsid.lower(), chromosome, position, genotype


def _native_load(path: Path, resolver: IidResolver | None) -> GenomeData | None:
    """Parse via the native Rust extension when it's installed, else ``None``.

    Falls back (returns ``None``) on any decode/IO error so the tolerant
    pure-Python reader can take over.
    """
    if _native is None:
        return None
    iid_map = resolver.as_dict() if resolver is not None else {}
    try:
        calls, total, diploid, internal, resolved = _native.parse_genome_file(
            str(path), iid_map
        )
    except OSError:
        return None
    if diploid == 0:
        raise ValueError(
            f"{path} contained no diploid SNP calls — is it a 23andMe or "
            "AncestryDNA raw file?"
        )
    return GenomeData(
        calls=calls,
        total_markers=total,
        diploid_markers=diploid,
        source=path,
        internal_ids=internal,
        resolved_ids=resolved,
    )


def load_genome(
    path: str | Path, *, resolver: IidResolver | None = None
) -> GenomeData:
    """Load and index a 23andMe or AncestryDNA genome file.

    Args:
        resolver: when given, 23andMe internal ``i`` ids are translated to their
            standard RSID so they can match SNPedia/GWAS entries.  The original
            ``i`` key is kept as an alias too, so lookups by either id succeed.

    Raises:
        FileNotFoundError: if ``path`` does not exist.
        ValueError: if the file has no usable diploid SNP calls (i.e. it is
            almost certainly not a 23andMe export).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Genome file not found: {path}")

    native = _native_load(path, resolver)
    if native is not None:
        return native

    calls: dict[str, str] = {}
    total = 0
    diploid = 0
    internal = 0
    resolved = 0
    for rsid, _chromosome, _position, genotype in _iter_data_lines(path):
        total += 1
        key = normalize_genotype(genotype)
        if key is None:
            continue
        diploid += 1
        marker = rsid.lower()
        calls.setdefault(marker, key)  # keep first; duplicates are rare
        if is_internal_id(marker):
            internal += 1
            canonical = resolver.resolve(marker) if resolver is not None else marker
            if canonical != marker:
                resolved += 1
                calls.setdefault(canonical, key)  # alias to the standard RSID

    if diploid == 0:
        raise ValueError(
            f"{path} contained no diploid SNP calls — is it a 23andMe or "
            "AncestryDNA raw file?"
        )

    return GenomeData(
        calls=calls,
        total_markers=total,
        diploid_markers=diploid,
        source=path,
        internal_ids=internal,
        resolved_ids=resolved,
    )
