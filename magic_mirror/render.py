"""Compose the genetic overlay onto a video frame using Pillow.

This module is deliberately camera-agnostic: it takes and returns plain RGB
numpy arrays, so it can decorate a live webcam frame *or* a still background
without importing OpenCV.  A single :class:`Renderer` caches fonts and the
vignette mask so the per-frame cost stays low.
"""

from __future__ import annotations

import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from . import theme
from .config import Settings
from .models import FaceBox, Trait


def _wrap(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont,
          max_width: int, max_lines: int) -> list[str]:
    """Greedy word-wrap that ellipsises once it runs out of lines."""
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = f"{current} {word}".strip()
        if draw.textlength(trial, font=font) <= max_width or not current:
            current = trial
        else:
            lines.append(current)
            current = word
            if len(lines) == max_lines:
                break
    if len(lines) < max_lines and current:
        lines.append(current)

    if lines:  # ellipsise the last line if we truncated the text
        consumed = sum(len(line.split()) for line in lines)
        if consumed < len(words):
            last = lines[-1]
            while draw.textlength(last + "…", font=font) > max_width and last:
                last = last[:-1].rstrip()
            lines[-1] = last + "…"
    return lines


def _ease(t: float) -> float:
    """Smoothstep easing for fades; clamps to [0, 1]."""
    t = max(0.0, min(1.0, t))
    return t * t * (3.0 - 2.0 * t)


class Renderer:
    """Stateful, reusable compositor for one display size / config."""

    def __init__(
        self,
        config: Settings,
        *,
        header: str,
        subheader: str,
        genes: dict[str, str] | None = None,
    ) -> None:
        self.config = config
        self.header = header
        self.subheader = subheader
        self.genes = genes or {}
        # Fonts scale a touch with card width so the design holds at any size.
        s = config.card_width
        self.f_title = theme.load_font("medium", 30)
        self.f_sub = theme.load_font("regular", 15)
        self.f_rsid = theme.load_font("mono", 13)
        self.f_body = theme.load_font("regular", max(15, s // 20))
        self.f_meta = theme.load_font("mono", 12)
        self._vignette_cache: tuple[tuple[int, int], Image.Image] | None = None

    # -- background mood -----------------------------------------------------
    def _mood(self, base: Image.Image) -> Image.Image:
        """Desaturate and vignette the reflection so overlay text reads."""
        cfg = self.config
        if cfg.desaturation > 0:
            gray = base.convert("L").convert("RGB")
            base = Image.blend(base, gray, cfg.desaturation)
        if cfg.vignette > 0:
            base = Image.alpha_composite(base.convert("RGBA"), self._vignette(base.size))
            base = base.convert("RGB")
        return base

    def _vignette(self, size: tuple[int, int]) -> Image.Image:
        if self._vignette_cache and self._vignette_cache[0] == size:
            return self._vignette_cache[1]
        w, h = size
        # Radial darkening toward the corners, built from a distance field.
        yy, xx = np.mgrid[0:h, 0:w]
        cx, cy = w / 2.0, h / 2.0
        dist = np.sqrt(((xx - cx) / cx) ** 2 + ((yy - cy) / cy) ** 2)
        strength = np.clip((dist - 0.55) / 0.85, 0.0, 1.0) ** 1.6
        alpha = (strength * 255 * self.config.vignette).astype(np.uint8)
        overlay = np.zeros((h, w, 4), dtype=np.uint8)
        overlay[..., 3] = alpha  # black with rising alpha
        image = Image.fromarray(overlay, "RGBA")
        self._vignette_cache = (size, image)
        return image

    # -- public API ----------------------------------------------------------
    def render(
        self,
        frame_rgb: np.ndarray,
        traits: list[Trait],
        *,
        time_s: float,
        face: FaceBox | None = None,
    ) -> np.ndarray:
        """Return ``frame_rgb`` with the genetic overlay composited on top."""
        base = Image.fromarray(frame_rgb).convert("RGB")
        base = self._mood(base).convert("RGBA")
        overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        self._draw_header(draw, base.size)
        if traits:
            self._draw_trait_pages(overlay, base.size, traits, time_s, face)

        composed = Image.alpha_composite(base, overlay).convert("RGB")
        return np.asarray(composed)

    # -- header --------------------------------------------------------------
    def _draw_header(self, draw: ImageDraw.ImageDraw, size: tuple[int, int]) -> None:
        m = self.config.margin
        title = " ".join(self.header.upper())  # letter-spaced caps
        draw.text((m, m - 4), title, font=self.f_title, fill=theme.INK)
        draw.text((m, m + 32), self.subheader, font=self.f_sub, fill=theme.INK_MUTED)

    # -- trait cards ---------------------------------------------------------
    def _draw_trait_pages(
        self,
        overlay: Image.Image,
        size: tuple[int, int],
        traits: list[Trait],
        time_s: float,
        face: FaceBox | None,
    ) -> None:
        cfg = self.config
        per_page = cfg.cards_per_page
        pages = max(1, math.ceil(len(traits) / per_page))

        page = int(time_s // cfg.page_seconds) % pages
        local_t = time_s - (time_s // cfg.page_seconds) * cfg.page_seconds
        fade = cfg.fade_seconds
        alpha_in = _ease(local_t / fade)
        alpha_out = _ease((cfg.page_seconds - local_t) / fade)
        page_alpha = min(alpha_in, alpha_out)
        if page_alpha <= 0.01:
            return

        w, h = size
        card_w = cfg.card_width
        cards = traits[page * per_page : (page + 1) * per_page]

        face_center = face.center if face else (w // 2, h // 2)
        for index, trait in enumerate(cards):
            col = index % cfg.columns
            row = index // cfg.columns
            x = cfg.margin if col == 0 else w - cfg.margin - card_w
            # Start cards below the header, stagger each slightly for rhythm.
            y = cfg.margin + 78 + row * (self._card_height() + cfg.card_gap)
            # Per-card stagger so a page assembles rather than snapping in.
            stagger = _ease((local_t - index * 0.06) / fade)
            card_alpha = min(page_alpha, stagger if alpha_in < 1 else page_alpha)
            self._draw_connector(overlay, (x, y), card_w, face_center, card_alpha)
            self._draw_card(overlay, trait, x, y, card_w, card_alpha)

    def _card_height(self) -> int:
        # rsid line + up to 3 body lines + meter, with padding.
        return 16 + 20 + 3 * (self.f_body.size + 6) + 22 + 20

    def _draw_connector(
        self,
        overlay: Image.Image,
        card_xy: tuple[int, int],
        card_w: int,
        face_center: tuple[int, int],
        alpha: float,
    ) -> None:
        x, y = card_xy
        # Anchor on the card edge nearest the face.
        anchor_x = x + card_w if x < face_center[0] else x
        anchor_y = y + self._card_height() // 2
        line = Image.new("RGBA", overlay.size, (0, 0, 0, 0))
        d = ImageDraw.Draw(line)
        r, g, b, a = theme.CONNECTOR
        d.line(
            [(anchor_x, anchor_y), face_center],
            fill=(r, g, b, int(a * alpha)),
            width=1,
        )
        d.ellipse(
            [anchor_x - 3, anchor_y - 3, anchor_x + 3, anchor_y + 3],
            fill=(r, g, b, int(min(255, a * 2) * alpha)),
        )
        overlay.alpha_composite(line)

    def _draw_card(
        self,
        overlay: Image.Image,
        trait: Trait,
        x: int,
        y: int,
        width: int,
        alpha: float,
    ) -> None:
        height = self._card_height()
        pad = 16
        accent = theme.accent_for(trait.category)

        tile = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        d = ImageDraw.Draw(tile)

        # Card surface + hairline border + left accent bar.
        d.rounded_rectangle(
            [0, 0, width - 1, height - 1],
            radius=14,
            fill=theme.CARD_FILL,
            outline=theme.CARD_BORDER,
            width=1,
        )
        d.rounded_rectangle([0, 0, 4, height - 1], radius=2, fill=(*accent, 255))

        # gene + rsid + genotype (mono) as an eyebrow line; the gene leads so the
        # card reads like "ACTN3 · rs1815739 · C;C".
        gene = self.genes.get(trait.rsid, "")
        eyebrow = (
            f"{gene}  ·  {trait.rsid}  ·  {trait.display_genotype}"
            if gene
            else f"{trait.rsid}  ·  {trait.display_genotype}"
        )
        d.text((pad, 12), eyebrow, font=self.f_rsid, fill=(*accent, 255))

        # Repute glyph (Promethease-style +/−/•) in the top-right corner.
        glyph = theme.REPUTE_GLYPHS.get(trait.repute, "•")
        rep_color = theme.repute_color(trait.repute)
        gw = d.textlength(glyph, font=self.f_title)
        d.text((width - pad - gw, 6), glyph, font=self.f_title, fill=(*rep_color, 255))

        # Summary body, wrapped to at most three lines.
        body_y = 12 + 22
        lines = _wrap(d, trait.summary, self.f_body, width - 2 * pad, 3)
        for line in lines:
            d.text((pad, body_y), line, font=self.f_body, fill=theme.INK)
            body_y += self.f_body.size + 6

        # Magnitude meter: five pips + numeric, bottom-aligned.
        self._draw_meter(d, trait.magnitude, pad, height - 26, accent)

        # Apply the fade by scaling the tile's alpha channel, then composite.
        if alpha < 0.999:
            faded = tile.getchannel("A").point(lambda p: int(p * alpha))
            tile.putalpha(faded)
        overlay.alpha_composite(tile, (x, y))

    def _draw_meter(
        self,
        d: ImageDraw.ImageDraw,
        magnitude: float,
        x: int,
        y: int,
        accent: tuple[int, int, int],
    ) -> None:
        pips = 5
        filled = max(0, min(pips, round(magnitude)))
        pip_w, pip_h, gap = 22, 4, 6
        for i in range(pips):
            px = x + i * (pip_w + gap)
            on = i < filled
            color = (*accent, 255) if on else (*theme.INK_MUTED, 70)
            d.rounded_rectangle([px, y, px + pip_w, y + pip_h], radius=2, fill=color)
        label_x = x + pips * (pip_w + gap) + 4
        d.text((label_x, y - 6), f"{magnitude:g}", font=self.f_meta, fill=theme.INK_MUTED)
