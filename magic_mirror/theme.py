"""Colours, fonts and small drawing helpers for the mirror's look.

The palette is tuned for legibility over a live, slightly-darkened camera feed:
soft desaturated accents on near-white text, no pure primaries.
"""

from __future__ import annotations

from functools import lru_cache

from PIL import ImageFont

RGB = tuple[int, int, int]
RGBA = tuple[int, int, int, int]

# --- text / surfaces ---
INK: RGB = (244, 246, 250)          # near-white primary text
INK_MUTED: RGB = (168, 176, 190)    # secondary text (rsid, meta)
CARD_FILL: RGBA = (14, 16, 22, 190)  # translucent dark card
CARD_BORDER: RGBA = (255, 255, 255, 28)
CONNECTOR: RGBA = (255, 255, 255, 40)

# --- category accents (desaturated, all readable on dark) ---
ACCENTS: dict[str, RGB] = {
    "longevity": (214, 178, 92),    # gold
    "cognition": (162, 148, 226),   # violet
    "performance": (232, 156, 92),  # amber
    "protection": (110, 198, 168),  # teal-green
    "risk": (226, 128, 122),        # coral
    "pharma": (120, 202, 190),      # mint (drug response)
    "appearance": (128, 178, 224),  # sky
    "metabolism": (198, 160, 210),  # orchid
    "trait": (150, 162, 182),       # slate (default)
}


def accent_for(category: str) -> RGB:
    return ACCENTS.get(category, ACCENTS["trait"])


# --- repute (Promethease-style Good / Bad / Neutral) ---
REPUTE_COLORS: dict[str, RGB] = {
    "Good": (110, 198, 168),     # green
    "Bad": (226, 128, 122),      # coral-red
    "Neutral": (150, 162, 182),  # slate
}
REPUTE_GLYPHS = {"Good": "+", "Bad": "−", "Neutral": "•"}


def repute_color(repute: str) -> RGB:
    return REPUTE_COLORS.get(repute, REPUTE_COLORS["Neutral"])


# Font candidates in preference order, spanning macOS, common Linux and a
# Pillow-bundled last resort so the mirror never crashes for lack of a font.
_REGULAR_CANDIDATES = (
    "/System/Library/Fonts/HelveticaNeue.ttc",
    "/System/Library/Fonts/Helvetica.ttc",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "DejaVuSans.ttf",  # Pillow ships this; resolved by name
)
_MEDIUM_CANDIDATES = (
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/HelveticaNeue.ttc",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "DejaVuSans-Bold.ttf",
)
_MONO_CANDIDATES = (
    "/System/Library/Fonts/Menlo.ttc",
    "/System/Library/Fonts/Supplemental/Courier New.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    "DejaVuSansMono.ttf",
)

_FAMILIES = {
    "regular": _REGULAR_CANDIDATES,
    "medium": _MEDIUM_CANDIDATES,
    "mono": _MONO_CANDIDATES,
}


@lru_cache(maxsize=64)
def load_font(weight: str, size: int) -> ImageFont.FreeTypeFont:
    """Return a cached TrueType font for ``weight`` ('regular'|'medium'|'mono').

    Falls back down the candidate list, and finally to Pillow's built-in bitmap
    font so rendering can never fail outright.
    """
    for candidate in _FAMILIES.get(weight, _REGULAR_CANDIDATES):
        try:
            return ImageFont.truetype(candidate, size)
        except (OSError, ValueError):
            continue
    return ImageFont.load_default()
