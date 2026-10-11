from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from magic_mirror.models import GenomeData, Trait
from magic_mirror.report import (
    render_html,
    report_counts,
    trait_to_row,
    write_csv,
)

_GENOME = GenomeData(calls={}, total_markers=1000, diploid_markers=900, source=Path("g.txt"))


def _trait(rsid, geno, mag, repute, category="risk", summary="a thing"):
    return Trait(rsid=rsid, genotype=geno, magnitude=mag, summary=summary,
                 category=category, repute=repute)


def test_trait_to_row_shape_and_links():
    row = trait_to_row(_trait("rs53576", "GG", 2.5, "Good", "cognition", "calm"))
    assert row["rsid"] == "rs53576"
    assert row["genotype"] == "G;G"
    assert row["repute"] == "Good"
    assert row["url"].endswith("/rs53576")
    assert "Rs53576(G;G)" in row["gurl"]   # genotype page is capitalised


def test_report_counts():
    traits = [
        _trait("rs1", "AA", 3, "Good"),
        _trait("rs2", "CC", 2, "Bad"),
        _trait("rs3", "TT", 1, "Bad"),
        _trait("rs4", "GG", 0, "Neutral"),
    ]
    assert report_counts(traits) == {"Good": 1, "Bad": 2, "Neutral": 1, "total": 4}


def test_render_html_is_self_contained_and_embeds_data():
    traits = [_trait("rs1", "AA", 5, "Bad", summary="higher risk of x")]
    html = render_html(traits, genome=_GENOME, title="Genetic Reflexions")

    # No unresolved template tokens, and no external resources.
    assert "__" not in re.sub(r"[A-Za-z]", "", html) or "__DATA__" not in html
    assert "__DATA__" not in html and "__TOTAL__" not in html
    assert "http-equiv" not in html
    assert "cdn" not in html.lower() and "<script src=" not in html

    # Embedded JSON round-trips and carries our row.
    match = re.search(r"const DATA = (\[.*?\]);\nconst list", html, re.S)
    data = json.loads(match.group(1).replace("<\\/", "</"))
    assert data[0]["rsid"] == "rs1" and data[0]["repute"] == "Bad"
    assert "1,000 markers read" in html


def test_render_html_sorts_by_magnitude_desc():
    traits = [_trait("rslo", "AA", 1, "Good"), _trait("rshi", "CC", 9, "Bad")]
    html = render_html(traits, genome=_GENOME)
    match = re.search(r"const DATA = (\[.*?\]);\nconst list", html, re.S)
    data = json.loads(match.group(1).replace("<\\/", "</"))
    assert [d["rsid"] for d in data] == ["rshi", "rslo"]


def test_write_csv_roundtrip(tmp_path):
    traits = [_trait("rs1", "AA", 5, "Bad", "risk", "higher risk")]
    out = write_csv(traits, tmp_path / "sub" / "report.csv", {"rs1": "MCM6"})
    rows = list(csv.reader(out.read_text().splitlines()))
    assert rows[0] == ["rsid", "gene", "genotype", "magnitude", "repute",
                       "category", "summary", "url"]
    assert rows[1][0] == "rs1" and rows[1][1] == "MCM6" and rows[1][4] == "Bad"


def test_render_html_neutralises_closing_script_in_payload():
    # A summary containing </script> must not be able to break out of the
    # embedded data <script> block.
    traits = [_trait("rs1", "AA", 2, "Neutral", summary="x</script><b>evil")]
    html = render_html(traits, genome=_GENOME)
    payload = re.search(r"const DATA = (\[.*?\]);\nconst list", html, re.S).group(1)
    assert "</script>" not in payload      # closing tag defused
    assert "<\\/script>" in payload        # ...to its escaped form
