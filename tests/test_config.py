from __future__ import annotations

import pytest
from pydantic import ValidationError

from magic_mirror.config import Settings, Source


def test_defaults():
    s = Settings()
    assert s.source is Source.LIVE
    assert s.min_magnitude == 2.0
    assert s.cards_per_page == s.columns * s.rows_per_column
    assert s.cache_ttl_seconds == s.cache_ttl_days * 86_400


def test_env_override(monkeypatch):
    monkeypatch.setenv("MIRROR_MIN_MAGNITUDE", "3.5")
    monkeypatch.setenv("MIRROR_SOURCE", "csv")
    monkeypatch.setenv("MIRROR_DEVICE_ID", "2")
    s = Settings()
    assert s.min_magnitude == 3.5
    assert s.source is Source.CSV
    assert s.device_id == 2


@pytest.mark.parametrize(
    "kwargs",
    [{"min_magnitude": -1}, {"max_traits": 0}, {"face_smoothing": 2.0}, {"vignette": -0.5}],
)
def test_validation_rejects_bad_values(kwargs):
    with pytest.raises(ValidationError):
        Settings(**kwargs)
