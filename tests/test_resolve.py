from __future__ import annotations

import pytest

from magic_mirror import snpedia
from magic_mirror.config import Settings, Source
from magic_mirror.snpedia.live import SnpediaUnavailable


def test_resolve_table_csv_source():
    settings = Settings(source=Source.CSV, genome_from_usb=False)
    table = snpedia.resolve_table(settings, log=lambda *_: None)
    assert table  # loaded from the bundled result.csv


def test_resolve_table_auto_falls_back_to_csv(tmp_path, monkeypatch):
    def boom(*_a, **_k):
        raise SnpediaUnavailable("offline")

    monkeypatch.setattr(snpedia, "fetch_live", boom)
    settings = Settings(source=Source.AUTO, cache_dir=tmp_path, genome_from_usb=False)
    table = snpedia.resolve_table(settings, log=lambda *_: None)
    assert table  # fell back to the CSV snapshot


def test_resolve_table_live_reraises_without_snapshot(tmp_path, monkeypatch):
    def boom(*_a, **_k):
        raise SnpediaUnavailable("offline")

    monkeypatch.setattr(snpedia, "fetch_live", boom)
    settings = Settings(
        source=Source.LIVE,
        cache_dir=tmp_path,
        snpedia_path=tmp_path / "missing.csv",  # no snapshot to fall back to
        genome_from_usb=False,
    )
    with pytest.raises(SnpediaUnavailable):
        snpedia.resolve_table(settings, log=lambda *_: None)
