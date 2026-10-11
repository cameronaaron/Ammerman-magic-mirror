"""Orchestration: resolve inputs, then run the live mirror or a still poster.

The poster and reflection-building paths need only Pillow + numpy; the live
mirror additionally needs OpenCV (imported lazily).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageOps

from . import snpedia
from .config import Settings
from .genome import load_genome
from .iids import IidResolver
from .models import FaceBox, GenomeData, Trait
from .render import Renderer
from .sources import find_local_genome, find_usb_genome
from .traits import build_traits


@dataclass(slots=True)
class Reflection:
    """Everything the display needs: ranked traits and their provenance."""

    traits: list[Trait]
    genome: GenomeData
    matched_total: int  # your genotypes present in the loaded SNPedia table

    @property
    def subheader(self) -> str:
        markers = f"{self.genome.total_markers:,}"
        return f"{markers} markers read  ·  {len(self.traits)} notable traits"


def resolve_genome_path(settings: Settings, *, log=print) -> Path:
    """Pick a genome: inserted USB first, then the local drop folder, else default."""
    if settings.genome_from_usb:
        usb = find_usb_genome()
        if usb is not None:
            log(f"Found genome on removable media: {usb}")
            return usb
    local = find_local_genome(settings.genome_inbox)
    if local is not None:
        log(f"Found genome in local drop folder: {local}")
        return local
    return settings.genome_path


def build_resolver(settings: Settings) -> IidResolver | None:
    """The i-id → RSID resolver for this run, or ``None`` when disabled."""
    if not settings.resolve_iids:
        return None
    return IidResolver.from_files(settings.iid_map_path)


def build_all_traits(
    settings: Settings, *, log=print
) -> tuple[GenomeData, list[Trait]]:
    """Load the genome and return *every* SNPedia match with a summary.

    Unlike :func:`build_reflection` (which keeps only the notable top slice for
    the mirror), this is the full Promethease-style set the HTML report needs.
    """
    genome = load_genome(
        resolve_genome_path(settings, log=log), resolver=build_resolver(settings)
    )
    table = snpedia.resolve_table(settings, log=log)
    traits = build_traits(genome, table, min_magnitude=0.0, require_summary=True)
    return genome, traits


def build_reflection(
    settings: Settings, *, positive_only: bool = False, log=print
) -> Reflection:
    """Load genome + SNPedia and rank the user's notable traits.

    ``positive_only`` keeps just the 'good news' traits — used for shareable
    graphics that should flatter, not alarm.
    """
    genome = load_genome(
        resolve_genome_path(settings, log=log), resolver=build_resolver(settings)
    )
    table = snpedia.resolve_table(settings, log=log)
    all_matched = build_traits(genome, table, min_magnitude=0.0, require_summary=True)
    notable = build_traits(
        genome,
        table,
        min_magnitude=settings.min_magnitude,
        require_summary=True,
        positive_only=positive_only,
        limit=settings.max_traits,
    )
    return Reflection(traits=notable, genome=genome, matched_total=len(all_matched))


# --------------------------------------------------------------------------- #
# Live mirror
# --------------------------------------------------------------------------- #
def run_live(settings: Settings, *, fullscreen: bool = False, log=print) -> None:
    """Open the webcam and run the magic mirror until Esc/q is pressed."""
    import cv2  # local import: only the live path depends on OpenCV

    from .camera import Camera
    from .detect import create_detector

    reflection = build_reflection(settings, log=log)
    renderer = Renderer(
        settings, header=settings.window_title, subheader=reflection.subheader
    )
    detector = create_detector(settings, log=log)

    window = settings.window_title
    cv2.namedWindow(window, cv2.WINDOW_NORMAL)
    if fullscreen:
        cv2.setWindowProperty(window, cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_FULLSCREEN)

    log(f"✦ {reflection.subheader}")
    log("  Live mirror running — Esc/q to quit, s to save a snapshot.")

    start = time.monotonic()
    with Camera(settings.device_id, mirror=settings.mirror) as camera:
        while True:
            frame_rgb = camera.read_rgb()
            if frame_rgb is None:
                if cv2.waitKey(15) & 0xFF in (27, ord("q")):
                    break
                continue

            face = detector.detect(frame_rgb[:, :, ::-1])  # detector wants BGR
            composed = renderer.render(
                frame_rgb, reflection.traits, time_s=time.monotonic() - start, face=face
            )
            cv2.imshow(window, composed[:, :, ::-1])  # back to BGR for OpenCV

            key = cv2.waitKey(15) & 0xFF
            if key in (27, ord("q")):
                break
            if key == ord("s"):
                log(f"  Saved snapshot → {_save_snapshot(composed)}")

    cv2.destroyAllWindows()


def _save_snapshot(rgb: np.ndarray) -> Path:
    out = Path("output") / "mirror_snapshot.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb).save(out)
    return out


# --------------------------------------------------------------------------- #
# Shareable export (turn a photo + genome into an Instagram-ready image)
# --------------------------------------------------------------------------- #
#: Named output sizes for the ``share`` export, in pixels (w, h).
SHARE_SIZES: dict[str, tuple[int, int]] = {
    "square": (1080, 1080),     # Instagram feed (1:1)
    "portrait": (1080, 1350),   # Instagram feed (4:5) — the tallest feed post
    "story": (1080, 1920),      # Instagram / TikTok story (9:16)
    "wide": (1280, 800),        # landscape, matches the live window
}
DEFAULT_SHARE_FORMAT = "portrait"


def _fit_cover(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    """Scale ``image`` to *cover* ``size`` then centre-crop — no distortion."""
    tw, th = size
    w, h = image.size
    if w == 0 or h == 0:
        return image.resize(size)
    scale = max(tw / w, th / h)
    resized = image.resize((max(1, round(w * scale)), max(1, round(h * scale))))
    nw, nh = resized.size
    left, top = (nw - tw) // 2, (nh - th) // 2
    return resized.crop((left, top, left + tw, top + th))


def detect_face_in_image(
    rgb: np.ndarray, settings: Settings, *, log=print
) -> FaceBox | None:
    """Best-effort face box for a still photo; ``None`` if OpenCV is unavailable."""
    try:
        import cv2  # noqa: F401  — presence check only

        from .detect import create_detector
    except ImportError:
        return None
    try:
        detector = create_detector(settings, log=log)
    except RuntimeError:
        return None
    return detector.detect(rgb[:, :, ::-1])  # detector wants BGR


def compose_share_image(
    settings: Settings,
    reflection: Reflection,
    photo: Image.Image,
    *,
    fmt: str = DEFAULT_SHARE_FORMAT,
    detect: bool = True,
    genes: dict[str, str] | None = None,
    log=print,
) -> Image.Image:
    """Composite the genetic overlay onto ``photo`` at an Instagram-ready size."""
    size = SHARE_SIZES.get(fmt, SHARE_SIZES[DEFAULT_SHARE_FORMAT])
    work = settings.model_copy()

    # Honour the phone's EXIF orientation before cropping, else selfies rotate.
    photo = ImageOps.exif_transpose(photo)
    base = _fit_cover(photo.convert("RGB"), size)
    base_rgb = np.asarray(base)

    renderer = Renderer(
        work, header=work.window_title, subheader=reflection.subheader, genes=genes
    )
    # Grow the number of card rows to fill taller formats (story shows more).
    available = size[1] - (work.margin + 78) - work.margin
    work.rows_per_column = max(1, available // (renderer._card_height() + work.card_gap))

    face = detect_face_in_image(base_rgb, work, log=log) if detect else None
    if face is None:
        # Upper-centre placeholder — where a portrait subject's head usually is.
        fw = size[0] // 4
        face = FaceBox(size[0] // 2 - fw // 2, size[1] // 4, fw, fw)

    composed = renderer.render(
        base_rgb, reflection.traits, time_s=work.page_seconds / 2, face=face
    )
    return Image.fromarray(composed)


def run_share(
    settings: Settings,
    photo: Path,
    output: Path,
    *,
    fmt: str = DEFAULT_SHARE_FORMAT,
    detect: bool = True,
    positive_only: bool = False,
    log=print,
) -> Path:
    """Render a shareable image from ``photo`` + the resolved genome to ``output``."""
    if fmt not in SHARE_SIZES:
        raise ValueError(
            f"unknown format '{fmt}'; choose one of {', '.join(SHARE_SIZES)}"
        )
    if not Path(photo).exists():
        raise FileNotFoundError(f"Photo not found: {photo}")
    reflection = build_reflection(settings, positive_only=positive_only, log=log)
    if positive_only:
        # Re-rank the good-news traits so the graphic leads with the coolest,
        # dropping near-duplicate summaries.
        from .traits import rank_cool

        reflection.traits = rank_cool(reflection.traits)

    # Label the cards with their gene (ACTN3, OXTR, …) — cached, best-effort.
    genes: dict[str, str] = {}
    try:
        from .snpedia.genes import fetch_genes

        genes = fetch_genes(
            [t.rsid for t in reflection.traits],
            cache_dir=settings.cache_dir,
            timeout=settings.http_timeout,
            log=log,
        )
    except Exception:  # noqa: BLE001 - genes are a nicety; never block the graphic
        genes = {}

    image = compose_share_image(
        settings, reflection, Image.open(photo),
        fmt=fmt, detect=detect, genes=genes, log=log,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output)
    return output


# --------------------------------------------------------------------------- #
# Still poster (no camera / OpenCV required)
# --------------------------------------------------------------------------- #
def studio_background(size: tuple[int, int], face: FaceBox) -> np.ndarray:
    """A quiet gradient with a soft glow standing in for a reflection."""
    w, h = size
    top = np.array([26, 28, 34], dtype=np.float64)
    bottom = np.array([12, 13, 17], dtype=np.float64)
    ramp = np.linspace(0, 1, h)[:, None, None]
    grad = (top * (1 - ramp) + bottom * ramp).astype(np.uint8)
    canvas = np.repeat(grad, w, axis=1)

    yy, xx = np.mgrid[0:h, 0:w]
    cx, cy = face.center
    dist = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2) / (0.6 * max(face.w, face.h))
    glow = np.clip(1 - dist, 0, 1)[..., None] ** 2
    return np.clip(canvas + glow * np.array([40, 44, 58]), 0, 255).astype(np.uint8)


def run_poster(
    settings: Settings,
    output: Path,
    *,
    size: tuple[int, int] = (1280, 800),
    at_seconds: float = 1.2,
    background: Path | None = None,
    log=print,
) -> Path:
    """Render a single decorated frame to ``output`` and return its path."""
    reflection = build_reflection(settings, log=log)
    renderer = Renderer(
        settings, header=settings.window_title, subheader=reflection.subheader
    )

    face = FaceBox(size[0] // 2 - 130, size[1] // 2 - 60, 260, 300)
    if background is not None:
        base = np.asarray(Image.open(background).convert("RGB").resize(size))
    else:
        base = studio_background(size, face)

    composed = renderer.render(base, reflection.traits, time_s=at_seconds, face=face)
    output.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(composed).save(output)
    return output
