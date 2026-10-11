from __future__ import annotations

import httpx
import respx

from magic_mirror.snpedia import genes
from magic_mirror.snpedia.genes import _extract_gene, fetch_genes


def test_extract_gene_from_wikitext():
    assert _extract_gene("{{Rsnum\n|rsid=53576\n|Gene=OXTR\n|Chromosome=3\n}}") == "OXTR"
    assert _extract_gene("| Gene = ACTN3 |") == "ACTN3"
    assert _extract_gene("no gene here") == ""


def _page(title: str, content: str) -> dict:
    return {"title": title, "revisions": [{"slots": {"main": {"content": content}}}]}


@respx.mock
def test_fetch_genes_parses_and_caches(tmp_path):
    respx.get(genes.API_ENDPOINT).mock(
        return_value=httpx.Response(
            200,
            json={"query": {"pages": [
                _page("Rs53576", "{{Rsnum|Gene=OXTR|Chromosome=3}}"),
                _page("Rs1815739", "{{Rsnum|Gene=ACTN3}}"),
                _page("Rs999", "no gene template"),
            ]}},
        )
    )
    result = fetch_genes(
        ["rs53576", "rs1815739", "rs999", "i123"],  # i-id ignored
        cache_dir=tmp_path,
        log=lambda *_: None,
    )
    assert result == {"rs53576": "OXTR", "rs1815739": "ACTN3"}  # rs999 has no gene
    assert (tmp_path / "snpedia_genes.json").exists()


@respx.mock
def test_fetch_genes_uses_cache_second_time(tmp_path):
    route = respx.get(genes.API_ENDPOINT).mock(
        return_value=httpx.Response(
            200, json={"query": {"pages": [_page("Rs1", "{{Rsnum|Gene=AAA}}")]}}
        )
    )
    kwargs = dict(cache_dir=tmp_path, log=lambda *_: None)
    fetch_genes(["rs1"], **kwargs)
    fetch_genes(["rs1"], **kwargs)  # served from cache, no second request
    assert route.call_count == 1


@respx.mock
def test_fetch_genes_network_failure_is_graceful(tmp_path):
    respx.get(genes.API_ENDPOINT).mock(return_value=httpx.Response(503))
    # Should not raise; returns empty (nothing cached yet).
    assert fetch_genes(["rs1"], cache_dir=tmp_path, log=lambda *_: None) == {}
