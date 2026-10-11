from __future__ import annotations

import httpx
import respx

from magic_mirror.snpedia.gwas import GwasAssoc, _parse_response, fetch_gwas


def _assoc(trait, allele, raf, pval):
    return {
        "pvalue": pval,
        "riskFrequency": raf,
        "orPerCopyNum": None,
        "betaNum": 0.05,
        "efoTraits": [{"trait": trait}],
        "loci": [{"strongestRiskAlleles": [{"riskAlleleName": f"rs1815739-{allele}"}]}],
    }


def test_parse_response_extracts_best_per_trait():
    data = {"_embedded": {"associations": [
        _assoc("muscle", "C", "0.55", 2e-25),
        _assoc("muscle", "C", "0.55", 8e-12),   # same trait/allele, weaker p
        _assoc("odd", "?", "NR", 1e-5),          # non-ACGT allele -> skipped
    ]}}
    rows = _parse_response("rs1815739", data)
    assert len(rows) == 1                        # one (trait, allele) kept
    trait, allele, raf, pval, _effect = rows[0]
    assert (trait, allele, raf) == ("muscle", "C", 0.55)
    assert pval == 2e-25                          # the stronger association


def test_gwas_assoc_dosage():
    a = GwasAssoc("rs1", "muscle", "C", 0.5, 1e-9, "0.05")
    assert a.dosage("CC") == 2
    assert a.dosage("CG") == 1
    assert a.dosage("GG") == 0


@respx.mock
def test_fetch_gwas_parses_and_caches(tmp_path):
    respx.get(url__regex=r".*rs1815739.*").mock(
        return_value=httpx.Response(200, json={"_embedded": {"associations": [
            _assoc("muscle", "C", "0.55", 2e-25)]}})
    )
    result = fetch_gwas(["rs1815739", "i9"], cache_dir=tmp_path, log=lambda *_: None)
    assert "rs1815739" in result and "i9" not in result
    assoc = result["rs1815739"][0]
    assert assoc.risk_allele == "C" and assoc.dosage("CC") == 2
    assert (tmp_path / "gwas_associations.json").exists()


@respx.mock
def test_fetch_gwas_network_failure_is_graceful(tmp_path):
    respx.get(url__regex=r".*").mock(return_value=httpx.Response(503))
    assert fetch_gwas(["rs1"], cache_dir=tmp_path, log=lambda *_: None) == {}
