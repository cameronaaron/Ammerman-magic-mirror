from __future__ import annotations

from magic_mirror.genome import load_genome
from magic_mirror.iids import (
    DEFAULT_IID_MAP,
    IidResolver,
    is_internal_id,
    is_rsid,
    load_iid_map,
)


def test_id_shape_predicates():
    assert is_internal_id("i713426")
    assert is_internal_id("I713426")      # case-insensitive
    assert not is_internal_id("rs123")
    assert is_rsid("rs2942")
    assert not is_rsid("i713426")
    assert not is_internal_id("interleukin")  # not i + digits


def test_load_iid_map_parses_and_skips_junk(tmp_path):
    path = tmp_path / "map.tsv"
    path.write_text(
        "# header comment\n"
        "\n"
        "i713426\trs2942\n"      # tab-separated
        "i1,rs99\n"              # comma-separated
        "i2\tnope\n"             # bad rsid -> skipped
        "garbage line\n"         # no second column -> skipped
        "I3\tRS7\n"              # upper-cased -> lowered
    )
    mapping = load_iid_map(path)
    assert mapping == {"i713426": "rs2942", "i1": "rs99", "i3": "rs7"}


def test_load_iid_map_missing_file_is_empty(tmp_path):
    assert load_iid_map(tmp_path / "absent.tsv") == {}


def test_resolver_translates_only_known_internal_ids():
    resolver = IidResolver({"i713426": "rs2942"})
    assert resolver.resolve("i713426") == "rs2942"
    assert resolver.resolve("I713426") == "rs2942"   # normalised input
    assert resolver.resolve("i999") == "i999"        # unknown -> unchanged
    assert resolver.resolve("rs5") == "rs5"          # already an rsid
    assert "i713426" in resolver
    assert len(resolver) == 1


def test_resolver_unresolved_lists_only_missing_internal_ids():
    resolver = IidResolver({"i1": "rs1"})
    missing = resolver.unresolved(["i1", "i2", "i3", "rs9"])
    assert missing == ["i2", "i3"]


def test_from_files_later_source_wins(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    a.write_text("i1\trs1\n")
    b.write_text("i1\trs2\ni9\trs9\n")
    resolver = IidResolver.from_files(a, b)
    assert resolver.resolve("i1") == "rs2"   # b overrides a
    assert resolver.resolve("i9") == "rs9"


def test_shipped_crosswalk_loads():
    # The bundled file must at least parse without error (may be empty).
    assert isinstance(load_iid_map(DEFAULT_IID_MAP), dict)


def test_load_genome_resolves_internal_ids(tmp_path):
    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "rs1\t1\t100\tAG\n"
        "i713426\t1\t200\tCC\n"    # internal id, resolvable
        "i999\t1\t300\tTT\n"       # internal id, unknown
    )
    resolver = IidResolver({"i713426": "rs2942"})
    genome = load_genome(genome_file, resolver=resolver)

    assert genome.internal_ids == 2
    assert genome.resolved_ids == 1
    # Both the original i-id and its RSID alias resolve to the same call.
    assert genome.calls["i713426"] == "CC"
    assert genome.calls["rs2942"] == "CC"
    assert "i999" in genome.calls


def test_load_genome_without_resolver_keeps_internal_ids(tmp_path):
    genome_file = tmp_path / "g.txt"
    genome_file.write_text(
        " rsid\tchromosome\tposition\tgenotype\n"
        "i713426\t1\t200\tCC\n"
    )
    genome = load_genome(genome_file)
    assert genome.internal_ids == 1
    assert genome.resolved_ids == 0
    assert genome.calls["i713426"] == "CC"
