"""Join a parsed genome against SNPedia into ranked, displayable traits.

Each :class:`~magic_mirror.models.Trait` is tagged with a coarse *category* (used
only to pick an accent colour) inferred from keywords in its summary.
"""

from __future__ import annotations

from .models import GenomeData, SnpediaEntry, Trait

# Category -> keywords.  First match wins, so order most-specific to least.
_CATEGORY_KEYWORDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("longevity", ("longevity", "lifespan", "aging", "ageing", "nonagenarian")),
    ("cognition", ("memory", "optimis", "empath", "stress", "intelligence", "iq")),
    ("performance", ("muscle", "sprint", "power", "athlet", "endurance", "sleeper")),
    ("pharma", ("drug", "plavix", "warfarin", "detoxif", "response to treatment")),
    (
        "protection",
        ("reduced", "resistance", "decreased", "lower risk", "protect", "half the"),
    ),
    (
        "risk",
        ("risk", "cancer", "disease", "asthma", "gout", "infarction", "sclerosis",
         "diabetes", "carrier"),
    ),
    ("appearance", ("earwax", "hair", "eye", "skin", "odour", "odor", "height",
                    "bald", "light-skinned", "myopia")),
    ("metabolism", ("caffeine", "alcohol", "lactose", "metaboli", "diet", "weight",
                    "cholesterol")),
)
_DEFAULT_CATEGORY = "trait"


def categorize(summary: str) -> str:
    text = summary.lower()
    for category, keywords in _CATEGORY_KEYWORDS:
        if any(keyword in text for keyword in keywords):
            return category
    return _DEFAULT_CATEGORY


# Category -> repute when SNPedia leaves Repute unset (offline / CSV snapshot).
_CATEGORY_REPUTE = {"protection": "Good", "risk": "Bad"}


# Categories that read as flattering on a shareable graphic even when SNPedia
# leaves Repute neutral (athleticism, longevity, cognition, protection, looks).
_POSITIVE_CATEGORIES = {
    "performance", "longevity", "cognition", "protection", "appearance",
}


def is_positive(trait) -> bool:
    """True for a 'good news' trait: never Bad, and either Good or a nice category."""
    if trait.repute == "Bad":
        return False
    return trait.repute == "Good" or trait.category in _POSITIVE_CATEGORIES


# --------------------------------------------------------------------------- #
# "Coolness" ranking for shareable graphics — flatter, don't bore.
# --------------------------------------------------------------------------- #
_COOL_CATEGORY_BOOST = {
    "performance": 3.0, "longevity": 3.0, "cognition": 2.5,
    "protection": 1.5, "appearance": 0.5,
}
_COOL_KEYWORDS = {
    # Superpower flexes — weighted high so they lead the graphic.
    "sprint": 4.0, "muscle": 3.5, "athlet": 3.5, "power": 2.5, "endurance": 3.0,
    "sleeper": 3.0, "faster": 2.5, "better perform": 3.0,
    "lifespan": 4.0, "longevity": 4.0, "nonagenarian": 3.5, "centenarian": 4.0,
    "intelligence": 4.0, "iq": 3.5, "memory": 3.5, "cognit": 2.5,
    "empath": 2.5, "optimis": 2.5,
    "resistance": 3.5, "resistant": 3.5, "immune": 2.5,
    "increased": 1.5, "advantage": 2.5, "elite": 3.0, "protect": 1.0,
    "reduced risk": 0.5,  # good news, but clinical/repetitive — don't crowd out flexes
}
_MEH_KEYWORDS = {
    "earwax": -4.0, "odour": -3.0, "odor": -3.0, "ancestry": -3.5,
    "light-skinned": -3.0, "hair": -2.0, "bald": -1.5, "myopia": -1.5,
    "carrier": -2.0, "?": -1.2,       # speculative summaries end in "?"
    "chronic": -0.8, "disease": -0.6,  # generic clinical phrasing (still positive)
}


def cool_score(trait) -> float:
    """A heuristic 'how cool does this read on Instagram' score."""
    score = trait.magnitude + _COOL_CATEGORY_BOOST.get(trait.category, 0.0)
    if trait.repute == "Good":
        score += 1.0
    text = trait.summary.lower()
    for keyword, weight in _COOL_KEYWORDS.items():
        if keyword in text:
            score += weight
    for keyword, weight in _MEH_KEYWORDS.items():
        if keyword in text:
            score += weight
    return score


def rank_cool(traits: list[Trait], *, limit: int | None = None) -> list[Trait]:
    """Dedupe near-identical summaries and rank by coolness for a graphic."""
    best: dict[str, Trait] = {}
    for trait in traits:
        key = " ".join(trait.summary.lower().split())  # normalise whitespace/case
        if key not in best or cool_score(trait) > cool_score(best[key]):
            best[key] = trait
    ranked = sorted(best.values(), key=lambda t: (-cool_score(t), -t.magnitude, t.rsid))
    return ranked[:limit] if limit is not None else ranked


def resolve_repute(raw: str, category: str) -> str:
    """Normalise SNPedia's Repute to Good/Bad/Neutral, falling back on category.

    SNPedia's "Repute" field is authoritative when present; when it is blank
    (the frozen CSV carries none), we infer a coarse Good/Bad from the summary
    category so the report still colour-codes protective vs risk genotypes.
    """
    value = (raw or "").strip().lower()
    if value == "good":
        return "Good"
    if value == "bad":
        return "Bad"
    return _CATEGORY_REPUTE.get(category, "Neutral")


def build_traits(
    genome: GenomeData,
    snpedia: dict[tuple[str, str], SnpediaEntry],
    *,
    min_magnitude: float = 0.0,
    require_summary: bool = True,
    positive_only: bool = False,
    limit: int | None = None,
) -> list[Trait]:
    """Return the user's traits, ranked most-notable first.

    Args:
        min_magnitude: drop anything below this SNPedia magnitude.
        require_summary: drop entries with no human-readable summary.
        positive_only: keep only 'good news' traits (see :func:`is_positive`).
        limit: keep at most this many after ranking.
    """
    traits: list[Trait] = []
    for rsid, genotype in genome.calls.items():
        entry = snpedia.get((rsid, genotype))
        if entry is None or entry.magnitude < min_magnitude:
            continue
        if require_summary and not entry.summary:
            continue
        category = categorize(entry.summary)
        trait = Trait(
            rsid=rsid,
            genotype=genotype,
            magnitude=entry.magnitude,
            summary=entry.summary,
            category=category,
            repute=resolve_repute(entry.repute, category),
        )
        if positive_only and not is_positive(trait):
            continue
        traits.append(trait)

    # Highest magnitude first; break ties by rsid so the rotation is stable.
    traits.sort(key=lambda t: (-t.magnitude, t.rsid))
    return traits[:limit] if limit is not None else traits
