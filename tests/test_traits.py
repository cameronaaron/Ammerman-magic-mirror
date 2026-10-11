from __future__ import annotations

from magic_mirror.genome import load_genome
from magic_mirror.models import SnpediaEntry, Trait
from magic_mirror.snpedia.parse import load_csv
from magic_mirror.traits import (
    build_traits,
    categorize,
    cool_score,
    is_positive,
    rank_cool,
    resolve_repute,
)


def _write(genome_file, csv_file):
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs10\t1\t1\tTA\n"      # reversed order vs snpedia (AT) — must still match
        "rs20\t1\t2\tGG\n"
        "rs30\t1\t3\tCC\n"
    )
    csv_file.write_text(
        ",Magnitude,Summary\n"
        'Rs10(A;T),3.0,"High interest trait"\n'
        'Rs20(G;G),1.0,"Low magnitude"\n'
        'Rs30(C;C),4.0,""\n'         # no summary -> dropped when required
    )


def test_build_traits_ranks_and_filters(tmp_path):
    genome_file, csv_file = tmp_path / "g.txt", tmp_path / "s.csv"
    _write(genome_file, csv_file)
    genome, snpedia = load_genome(genome_file), load_csv(csv_file)

    traits = build_traits(genome, snpedia, min_magnitude=2.0, require_summary=True)
    assert [t.rsid for t in traits] == ["rs10"]
    assert traits[0].display_genotype == "A;T"   # order-independent match

    ranked = build_traits(genome, snpedia, min_magnitude=0.0)
    assert [t.rsid for t in ranked] == ["rs10", "rs20"]

    assert build_traits(genome, snpedia, min_magnitude=2.0, limit=0) == []


def test_categorize():
    assert categorize("Reduced risk of Breast Cancer") == "protection"
    assert categorize("Optimistic and empathetic") == "cognition"
    assert categorize("blue eye color, 99% of the time") == "appearance"
    assert categorize("something totally unrelated") == "trait"


def test_resolve_repute_prefers_snpedia_then_falls_back_to_category():
    assert resolve_repute("Good", "risk") == "Good"       # SNPedia wins
    assert resolve_repute("bad", "protection") == "Bad"   # case-insensitive
    assert resolve_repute("", "protection") == "Good"     # inferred from category
    assert resolve_repute("", "risk") == "Bad"
    assert resolve_repute("", "appearance") == "Neutral"  # no signal


def _trait(repute, category):
    return Trait("rs1", "AA", 2.0, "summary", category, repute)


def test_is_positive():
    assert is_positive(_trait("Good", "risk"))            # Good always positive
    assert not is_positive(_trait("Bad", "protection"))   # Bad never positive
    assert is_positive(_trait("Neutral", "performance"))  # flattering category
    assert not is_positive(_trait("Neutral", "risk"))     # neutral + risk => no


def test_cool_score_prefers_flattering_over_mundane():
    sprinter = Trait("rs1", "AA", 2.0, "Likely sprinter, better muscle", "performance", "Good")
    earwax = Trait("rs2", "AA", 2.0, "wet earwax", "appearance", "Neutral")
    assert cool_score(sprinter) > cool_score(earwax)


def test_rank_cool_dedupes_identical_summaries():
    traits = [
        Trait("rs1", "AA", 2.0, "associated with higher HDL cholesterol", "metabolism", "Neutral"),
        Trait("rs2", "AA", 2.0, "associated with higher HDL cholesterol", "metabolism", "Neutral"),
        Trait("rs3", "AA", 2.2, "Likely sprinter", "performance", "Good"),
    ]
    ranked = rank_cool(traits)
    summaries = [t.summary for t in ranked]
    assert summaries.count("associated with higher HDL cholesterol") == 1
    assert "Likely sprinter" in summaries


def test_build_traits_positive_only(tmp_path):
    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t1\tAA\n"      # Bad -> excluded
        "rs2\t1\t2\tGG\n"      # Good -> kept
        "rs3\t1\t3\tCC\n"      # neutral risk -> excluded
    )
    genome = load_genome(genome_file)
    table = {
        ("rs1", "AA"): SnpediaEntry("rs1", "AA", 3.0, "higher risk", repute="Bad"),
        ("rs2", "GG"): SnpediaEntry("rs2", "GG", 2.5, "protective", repute="Good"),
        ("rs3", "CC"): SnpediaEntry("rs3", "CC", 2.0, "increased cancer risk"),
    }
    kept = build_traits(genome, table, min_magnitude=0.0, positive_only=True)
    assert [t.rsid for t in kept] == ["rs2"]


def test_build_traits_threads_repute(tmp_path):
    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t1\tAA\n"
        "rs2\t1\t2\tGG\n"
    )
    genome = load_genome(genome_file)
    table = {
        ("rs1", "AA"): SnpediaEntry("rs1", "AA", 3.0, "higher risk", repute="Bad"),
        ("rs2", "GG"): SnpediaEntry("rs2", "GG", 2.0, "some appearance trait"),
    }
    traits = {t.rsid: t for t in build_traits(genome, table, min_magnitude=0.0)}
    assert traits["rs1"].repute == "Bad"       # from SNPedia
    assert traits["rs2"].repute == "Neutral"   # no repute, neutral category
