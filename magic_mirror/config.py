"""Validated, environment-aware configuration (pydantic-settings).

Every knob can be overridden by an ``MIRROR_``-prefixed environment variable,
e.g. ``MIRROR_MIN_MAGNITUDE=2.5`` or ``MIRROR_DEVICE_ID=1``.  The rest of the
codebase carries no magic numbers — they all live here.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_GENOME = PROJECT_ROOT / "sample" / "genome.txt"
DEFAULT_SNPEDIA = PROJECT_ROOT / "result.csv"
DEFAULT_IID_MAP = PROJECT_ROOT / "resources" / "iid_to_rsid.tsv"
CASCADE_DIR = PROJECT_ROOT / "cascadeFiles"


class Source(str, Enum):
    """Where trait data comes from."""

    LIVE = "live"  # SNPedia API (disk-cached), CSV fallback if offline
    CSV = "csv"    # frozen result.csv snapshot only
    AUTO = "auto"  # live, silently fall back to CSV on any failure


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="MIRROR_", env_file=".env", extra="ignore"
    )

    # --- data ---
    genome_path: Path = DEFAULT_GENOME
    #: When set, prefer a genome auto-discovered on inserted removable media.
    genome_from_usb: bool = True
    #: Local project folder to check for a 23andMe file (after USB, before the
    #: configured default) — drop a genome in here instead of using a USB stick.
    genome_inbox: Path = PROJECT_ROOT / "dna"
    source: Source = Source.LIVE
    snpedia_path: Path = DEFAULT_SNPEDIA  # CSV snapshot / offline fallback
    #: Translate 23andMe internal ``i`` ids to standard RSIDs on load.
    resolve_iids: bool = True
    iid_map_path: Path = DEFAULT_IID_MAP  # i -> rs crosswalk (TSV)
    min_magnitude: float = Field(2.0, ge=0)
    max_traits: int = Field(60, gt=0)
    #: Live fetches pull every genotype above this floor once and cache it; the
    #: display then filters locally by ``min_magnitude``.  Decoupling the two
    #: means changing ``min_magnitude`` reuses the same cache instead of
    #: re-scanning SNPedia.
    snpedia_fetch_floor: float = Field(2.0, ge=0)

    # --- live SNPedia fetch ---
    cache_dir: Path = PROJECT_ROOT / ".cache"
    cache_ttl_days: float = Field(30.0, ge=0)
    http_timeout: float = Field(60.0, gt=0)
    force_refresh: bool = False

    # --- camera / face detection ---
    device_id: int = 0
    face_cascade: Path = CASCADE_DIR / "haarcascade_frontalface_default.xml"
    #: EMA factor for smoothing the tracked face box (0 = frozen, 1 = raw).
    face_smoothing: float = Field(0.35, ge=0, le=1)
    mirror: bool = True  # flip horizontally so it behaves like a real mirror

    # --- layout ---
    window_title: str = "Genetic Reflexions"
    columns: int = Field(2, gt=0)
    rows_per_column: int = Field(3, gt=0)
    card_width: int = 340
    margin: int = 36
    card_gap: int = 18

    # --- animation ---
    page_seconds: float = Field(6.0, gt=0)
    fade_seconds: float = Field(0.9, gt=0)

    # --- mood ---
    desaturation: float = Field(0.35, ge=0, le=1)
    vignette: float = Field(0.55, ge=0, le=1)

    # --- web (optional upload/accounts server) ---
    web_data_dir: Path = PROJECT_ROOT / "web_data"
    web_host: str = "127.0.0.1"
    web_port: int = Field(8000, gt=0)

    @field_validator(
        "genome_path", "snpedia_path", "cache_dir", "face_cascade", "iid_map_path",
        "web_data_dir", "genome_inbox",
    )
    @classmethod
    def _expand(cls, value: Path) -> Path:
        return value.expanduser()

    @property
    def cards_per_page(self) -> int:
        return max(1, self.columns * self.rows_per_column)

    @property
    def cache_ttl_seconds(self) -> float:
        return self.cache_ttl_days * 86_400
