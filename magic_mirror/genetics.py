"""Pure genetics helpers shared by the genome and SNPedia layers.

Kept dependency-free and side-effect-free so every other module can import it
without cycles.
"""

from __future__ import annotations

import re

# A valid diploid SNP call is two bases from this alphabet.  No-calls ("--"),
# indels ("DI") and single-allele sex-chromosome calls are intentionally excluded.
BASES = frozenset("ACGT")

# A SNPedia genotype page title: an rsid (rs…) or internal id (i…) followed by a
# parenthesised diploid call, e.g. "Rs17822931(T;T)" or "I3002432(A;A)".
_TITLE_RE = re.compile(r"((?:rs|i)\d+)\(([ACGT]);([ACGT])\)", re.IGNORECASE)


def normalize_genotype(genotype: str) -> str | None:
    """Return a strand/order-independent key for a diploid call, else ``None``.

    23andMe reports alleles in reference-strand order (``"AG"``) while SNPedia may
    list either order; sorting the alleles lets the two join regardless.  ``"AG"``
    and ``"GA"`` both map to ``"AG"``.
    """
    g = genotype.strip().upper()
    if len(g) == 2 and g[0] in BASES and g[1] in BASES:
        return "".join(sorted(g))
    return None


def parse_genotype_title(title: str) -> tuple[str, str] | None:
    """Parse a SNPedia genotype title into ``(rsid, normalized_genotype)``.

    Returns ``None`` for indels, single-allele or otherwise malformed titles.
    """
    match = _TITLE_RE.match(title.strip())
    if not match:
        return None
    genotype = normalize_genotype(match.group(2) + match.group(3))
    if genotype is None:
        return None
    return match.group(1).lower(), genotype


def coerce_magnitude(value: object) -> float:
    """Coerce a magnitude cell (str, int, float, list, or None) to a float.

    Semantic MediaWiki returns each property as a list, so unwrap that too.
    """
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0
