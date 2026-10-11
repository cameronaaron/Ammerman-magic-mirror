from __future__ import annotations

import httpx
import pytest
import respx

from magic_mirror.snpedia import live
from magic_mirror.snpedia.live import SnpediaUnavailable, fetch_live


def _page(results: dict, cont: int | None) -> httpx.Response:
    body: dict = {"query": {"results": results}}
    if cont is not None:
        body["query-continue-offset"] = cont
    return httpx.Response(200, json=body)


def _geno(mag, summary, repute=None):
    printouts = {"Magnitude": [mag], "Summary": [summary]}
    if repute is not None:
        printouts["Repute"] = [repute]
    return {"printouts": printouts}


@pytest.fixture(autouse=True)
def _fast(monkeypatch):
    # Don't actually sleep between pages / retries during tests.
    monkeypatch.setattr(live, "_RETRY_BACKOFF", 0)
    monkeypatch.setattr(live, "_INTER_PAGE_PAUSE", 0)


@respx.mock
def test_fetch_live_paginates(tmp_path):
    respx.get(live.API_ENDPOINT).mock(
        side_effect=[
            _page({"Rs53576(G;G)": _geno(2.5, "Optimistic")}, cont=500),
            _page({"Rs1815739(C;C)": _geno(2.2, "Sprinter")}, cont=None),
        ]
    )
    table = fetch_live(2.0, cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)
    assert set(table) == {("rs53576", "GG"), ("rs1815739", "CC")}
    assert table[("rs53576", "GG")].summary == "Optimistic"


@respx.mock
def test_fetch_live_parses_repute(tmp_path):
    respx.get(live.API_ENDPOINT).mock(
        side_effect=[_page({"Rs1801133(A;A)": _geno(3, "MTHFR", "Bad")}, cont=None)]
    )
    table = fetch_live(2.0, cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)
    assert table[("rs1801133", "AA")].repute == "Bad"


def test_rows_to_table_tolerates_legacy_and_new_cache():
    # Pre-Repute caches wrote 4-field rows; new ones write 5.
    table = live._rows_to_table([
        ["rs1", "AA", 2.0, "old row"],                 # legacy: no repute
        ["rs2", "GG", 3.0, "new row", "Good"],         # new: with repute
    ])
    assert table[("rs1", "AA")].repute == ""
    assert table[("rs2", "GG")].repute == "Good"


@respx.mock
def test_fetch_live_keeps_partial_on_tail_failure(tmp_path):
    respx.get(live.API_ENDPOINT).mock(
        side_effect=[
            _page({"Rs53576(G;G)": _geno(2.5, "Optimistic")}, cont=500),
            httpx.Response(502),  # retried, then loop keeps the first page
            httpx.Response(502),
            httpx.Response(502),
        ]
    )
    table = fetch_live(2.0, cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)
    assert set(table) == {("rs53576", "GG")}


@respx.mock
def test_fetch_live_raises_when_totally_unreachable(tmp_path):
    respx.get(live.API_ENDPOINT).mock(return_value=httpx.Response(503))
    with pytest.raises(SnpediaUnavailable):
        fetch_live(2.0, cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)


def test_fetch_live_uses_cache(tmp_path, monkeypatch):
    calls = {"n": 0}

    def fake_download(min_magnitude, timeout, log):
        from magic_mirror.models import SnpediaEntry

        calls["n"] += 1
        return [SnpediaEntry("rs53576", "GG", 2.5, "Optimistic")]

    monkeypatch.setattr(live, "_download", fake_download)
    kwargs = dict(cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)
    first = fetch_live(2.0, **kwargs)
    second = fetch_live(2.0, **kwargs)  # served from disk cache
    assert first == second
    assert calls["n"] == 1  # network hit once, then cached


def test_fetch_live_refresh_bypasses_cache(tmp_path, monkeypatch):
    calls = {"n": 0}

    def fake_download(min_magnitude, timeout, log):
        from magic_mirror.models import SnpediaEntry

        calls["n"] += 1
        return [SnpediaEntry("rs1", "AA", 3.0, "x")]

    monkeypatch.setattr(live, "_download", fake_download)
    kwargs = dict(cache_dir=tmp_path, ttl_seconds=9999, log=lambda *_: None)
    fetch_live(2.0, **kwargs)
    fetch_live(2.0, force_refresh=True, **kwargs)
    assert calls["n"] == 2
