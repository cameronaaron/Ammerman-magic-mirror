"""Resolve 23andMe *internal* marker IDs (``i…``) to standard dbSNP RSIDs.

23andMe's custom Illumina chip reports a chunk of its markers under proprietary
``i``-prefixed identifiers (``i713426``) rather than the community-standard
Reference SNP cluster IDs (``rs713426``).  SNPedia — and every downstream trait
database — is keyed on RSIDs, so an ``i``-only marker silently fails to match
even when the underlying SNP is richly annotated.

The catch is that **the full ``i`` → ``rs`` translation is not publicly known**:
23andMe has never published a complete crosswalk, so this ships a *partial,
extensible* lookup rather than pretending to be authoritative.  Mappings live in
a plain TSV resource (:data:`DEFAULT_IID_MAP`) so the table can grow without a
code change, and callers can layer in their own file.  Whatever we can't resolve
is left untouched — never guessed at.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from pathlib import Path

#: Ships alongside the package; format documented in the file's own header.
DEFAULT_IID_MAP = Path(__file__).resolve().parent.parent / "resources" / "iid_to_rsid.tsv"

_IID_RE = re.compile(r"^i\d+$", re.IGNORECASE)
_RSID_RE = re.compile(r"^rs\d+$", re.IGNORECASE)


def is_internal_id(marker: str) -> bool:
    """True if ``marker`` is a 23andMe internal id (``i`` followed by digits)."""
    return bool(_IID_RE.match(marker.strip()))


def is_rsid(marker: str) -> bool:
    """True if ``marker`` is a standard dbSNP RSID (``rs`` followed by digits)."""
    return bool(_RSID_RE.match(marker.strip()))


def load_iid_map(path: str | Path = DEFAULT_IID_MAP) -> dict[str, str]:
    """Load an ``i``-id → RSID crosswalk from a TSV file.

    The file is ``iid<TAB>rsid`` per line; ``#`` comments and blank lines are
    ignored, keys/values are lower-cased, and rows that aren't a clean
    ``i…`` → ``rs…`` pair are skipped rather than trusted.  A missing file yields
    an empty map (the resolver then simply passes markers through).
    """
    path = Path(path)
    if not path.exists():
        return {}

    mapping: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"[\t,]", line)
            if len(parts) < 2:
                continue
            iid, rsid = parts[0].strip().lower(), parts[1].strip().lower()
            if is_internal_id(iid) and is_rsid(rsid):
                mapping.setdefault(iid, rsid)
    return mapping


class IidResolver:
    """Translate internal ids to RSIDs from one or more crosswalk sources.

    Built from the shipped map by default; extend it with a user-supplied file
    or an in-memory mapping.  Later sources win on conflict, so a caller's own
    table can override the bundled one.
    """

    def __init__(self, mapping: Mapping[str, str] | None = None) -> None:
        self._map: dict[str, str] = dict(mapping or {})

    @classmethod
    def from_files(cls, *paths: str | Path) -> IidResolver:
        resolver = cls()
        for path in paths:
            resolver.update(load_iid_map(path))
        return resolver

    @classmethod
    def default(cls) -> IidResolver:
        """The resolver backed only by the package's shipped crosswalk."""
        return cls(load_iid_map())

    def update(self, mapping: Mapping[str, str]) -> IidResolver:
        for iid, rsid in mapping.items():
            iid, rsid = iid.strip().lower(), rsid.strip().lower()
            if is_internal_id(iid) and is_rsid(rsid):
                self._map[iid] = rsid
        return self

    def as_dict(self) -> dict[str, str]:
        """A copy of the ``i``-id → RSID map (e.g. to hand to the native parser)."""
        return dict(self._map)

    def __len__(self) -> int:
        return len(self._map)

    def __contains__(self, marker: str) -> bool:
        return marker.strip().lower() in self._map

    def resolve(self, marker: str) -> str:
        """Return the RSID for ``marker`` if known, else ``marker`` unchanged.

        Only internal ids are ever rewritten; anything already an RSID (or an
        unknown ``i``-id) is returned as-is and lower-cased for consistency.
        """
        key = marker.strip().lower()
        return self._map.get(key, key)

    def unresolved(self, markers: Iterable[str]) -> list[str]:
        """The internal ids in ``markers`` we have no translation for (sorted)."""
        missing = {
            m.strip().lower()
            for m in markers
            if is_internal_id(m) and m.strip().lower() not in self._map
        }
        return sorted(missing)
