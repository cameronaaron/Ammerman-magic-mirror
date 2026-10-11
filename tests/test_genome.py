from __future__ import annotations

import pytest

from magic_mirror.genome import load_genome


def test_load_genome_counts(tmp_path):
    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        "# comment\n"
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t100\tAG\n"
        "rs2\t1\t200\t--\n"      # no-call: counted, not diploid
        "rs3\tX\t300\tGG\n"
    )
    genome = load_genome(genome_file)
    assert genome.total_markers == 3
    assert genome.diploid_markers == 2
    assert len(genome) == 2
    assert genome.calls["rs1"] == "AG"
    assert "rs2" not in genome.calls


def test_load_genome_ancestrydna_format(tmp_path):
    genome_file = tmp_path / "AncestryDNA.txt"
    genome_file.write_text(
        "#AncestryDNA raw data download\n"
        "rsid\tchromosome\tposition\tallele1\tallele2\n"
        "rs1\t1\t100\tA\tG\n"        # diploid, split alleles -> "AG"
        "rs2\t1\t200\t0\t0\n"        # AncestryDNA no-call
        "rs3\t23\t300\tC\tC\n"       # numeric chromosome code (X)
    )
    genome = load_genome(genome_file)
    assert genome.total_markers == 3
    assert genome.diploid_markers == 2
    assert genome.calls["rs1"] == "AG"   # alleles joined + sorted
    assert "rs2" not in genome.calls     # 0/0 dropped as a no-call
    assert genome.calls["rs3"] == "CC"


def test_native_parser_matches_pure_python(tmp_path, monkeypatch):
    """When the Rust extension is built, it must agree with the Python parser."""
    from magic_mirror import genome as gmod

    if gmod._native is None:
        pytest.skip("native genome_core extension not built")

    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t100\tAG\n"
        "rs2\t1\t200\t--\n"
        "i5\t1\t300\tCC\n"
        "rs3\tX\t400\tGG\n"
    )
    native = gmod.load_genome(genome_file)
    monkeypatch.setattr(gmod, "_native", None)   # force the pure-Python path
    pure = gmod.load_genome(genome_file)

    assert native.calls == pure.calls
    assert (native.total_markers, native.diploid_markers, native.internal_ids) == (
        pure.total_markers, pure.diploid_markers, pure.internal_ids
    )


def test_load_genome_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_genome(tmp_path / "nope.txt")


def test_load_genome_rejects_non_genome(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("just some prose, not a genome at all\n")
    with pytest.raises(ValueError):
        load_genome(bad)
