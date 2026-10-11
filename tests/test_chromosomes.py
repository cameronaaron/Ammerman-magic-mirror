from __future__ import annotations

from magic_mirror.chromosomes import (
    compare_genomes,
    summarize,
    write_discordant_csv,
    write_summary_csv,
)
from magic_mirror.iids import IidResolver
from magic_mirror.models import SnpediaEntry

_HEADER = " rsid\tchromosome\tposition\tgenotype\n"


def _write(path, rows: str):
    path.write_text(_HEADER + rows)
    return path


def test_summarize_orders_and_counts(tmp_path):
    genome = _write(
        tmp_path / "g.txt",
        "rs1\t2\t100\tAG\n"      # chr 2
        "rs2\t1\t200\tCC\n"      # chr 1
        "rs3\t1\t300\t--\n"      # chr 1, no-call
        "rsx\tX\t400\tTT\n"      # chr X
        "rsm\tMT\t500\tGG\n",    # chr MT
    )
    stats = summarize(genome)
    order = [s.chromosome for s in stats]
    assert order == ["1", "2", "X", "MT"]        # autosomes then X then MT

    by = {s.chromosome: s for s in stats}
    assert by["1"].markers == 2 and by["1"].diploid == 1   # one no-call dropped
    assert by["2"].markers == 1 and by["2"].diploid == 1


def test_summarize_matches_against_table_with_resolver(tmp_path):
    genome = _write(
        tmp_path / "g.txt",
        "rs1\t1\t100\tAG\n"
        "i77\t1\t200\tCC\n"       # internal id resolvable to rs900
        "rs2\t1\t300\tTT\n",      # not in table
    )
    table = {
        ("rs1", "AG"): SnpediaEntry("rs1", "AG", 3.0, "trait one"),
        ("rs900", "CC"): SnpediaEntry("rs900", "CC", 2.0, "trait two"),
    }
    resolver = IidResolver({"i77": "rs900"})
    stats = summarize(genome, table, resolver=resolver)
    assert stats[0].chromosome == "1"
    assert stats[0].matched == 2   # rs1 directly, i77 via resolver -> rs900


def test_summarize_ancestrydna_maps_numeric_chromosomes(tmp_path):
    genome = tmp_path / "AncestryDNA.txt"
    genome.write_text(
        "rsid\tchromosome\tposition\tallele1\tallele2\n"
        "rs1\t1\t100\tA\tG\n"       # autosome
        "rs2\t23\t200\tC\tC\n"      # X
        "rs3\t24\t300\tA\tA\n"      # Y
        "rs4\t25\t350\tT\tT\n"      # XY / pseudoautosomal
        "rs5\t26\t400\tT\tT\n"      # MT
    )
    stats = summarize(genome)
    assert [s.chromosome for s in stats] == ["1", "X", "Y", "XY", "MT"]
    assert all(s.diploid == 1 for s in stats)


def test_compare_genomes_buckets_shared_snps(tmp_path):
    a = _write(
        tmp_path / "a.txt",
        "rs1\t1\t100\tAG\n"       # shared, concordant
        "rs2\t1\t200\tCC\n"       # shared, discordant (order-independent? no: CC vs GG)
        "rs3\t2\t300\tTT\n"       # only in A
        "rs4\tX\t400\tGA\n",      # shared, concordant (GA == AG)
    )
    b = _write(
        tmp_path / "b.txt",
        "rs1\t1\t100\tGA\n"       # AG == GA -> concordant
        "rs2\t1\t200\tGG\n"       # CC != GG -> discordant
        "rs4\tX\t400\tAG\n"       # concordant
        "rs9\t9\t900\tTT\n",      # only in B
    )
    result = compare_genomes(a, b)
    assert result.shared == 3           # rs1, rs2, rs4
    assert result.concordant == 2       # rs1, rs4
    assert result.discordant == 1       # rs2
    assert [d.rsid for d in result.discordant_snps] == ["rs2"]
    d = result.discordant_snps[0]
    assert d.chromosome == "1" and d.genotype_a == "CC" and d.genotype_b == "GG"
    assert result.by_chromosome["1"] == (2, 1, 1)   # shared, conc, disc
    assert abs(result.concordance_rate - 2 / 3) < 1e-9


def test_write_csv_roundtrip(tmp_path):
    genome = _write(tmp_path / "g.txt", "rs1\t1\t100\tAG\n")
    stats = summarize(genome)
    out = write_summary_csv(stats, tmp_path / "sub" / "summary.csv")
    assert out.exists()
    assert "chromosome,markers,diploid,matched" in out.read_text().splitlines()[0]

    other = _write(tmp_path / "o.txt", "rs1\t1\t100\tGG\n")
    result = compare_genomes(genome, other)
    disc = write_discordant_csv(result, tmp_path / "disc.csv")
    lines = disc.read_text().splitlines()
    assert lines[0] == "rsid,chromosome,position,genotype_a,genotype_b"
    assert lines[1].startswith("rs1,1,100,")
