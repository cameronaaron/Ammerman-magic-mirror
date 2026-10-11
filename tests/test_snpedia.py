from __future__ import annotations

import pytest

from magic_mirror.snpedia.parse import load_csv


def test_load_csv_parses_key(tmp_path):
    csv_file = tmp_path / "s.csv"
    csv_file.write_text(
        ",Magnitude,Summary\n"
        'Rs17822931(T;T),2.5,"Dry earwax."\n'
        "Rs1(DEL),0,skip me\n"      # malformed key -> skipped
        "I3002432(A;A),6.7,internal id kept\n"
    )
    table = load_csv(csv_file)
    assert ("rs17822931", "TT") in table
    assert table[("rs17822931", "TT")].magnitude == 2.5
    assert ("i3002432", "AA") in table   # internal ids parse too
    assert len(table) == 2


def test_load_csv_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_csv(tmp_path / "nope.csv")


def test_load_csv_empty_raises(tmp_path):
    csv_file = tmp_path / "s.csv"
    csv_file.write_text(",Magnitude,Summary\nnonsense,,\n")
    with pytest.raises(ValueError):
        load_csv(csv_file)


def test_snapshot_round_trips_through_its_own_writer(tmp_path):
    """The writer puts Repute before Summary; the loader once read the
    summary by position and showed "Good"/"Bad" as every trait's text."""
    from magic_mirror.models import SnpediaEntry
    from magic_mirror.snpedia.parse import load_csv, write_snapshot_csv

    entry = SnpediaEntry("rs17822931", "TT", 2.5, "Dry earwax, no odour", "Good")
    table = load_csv(write_snapshot_csv([entry], tmp_path / "result.csv"))
    loaded = table[("rs17822931", "TT")]
    assert (loaded.summary, loaded.repute) == ("Dry earwax, no odour", "Good")


def test_snapshot_without_a_repute_column_still_loads(tmp_path):
    from magic_mirror.snpedia.parse import load_csv

    path = tmp_path / "old.csv"
    path.write_text(',Magnitude,Summary\nRs17822931(T;T),2.5,"Dry earwax"\n')
    loaded = load_csv(path)[("rs17822931", "TT")]
    assert (loaded.summary, loaded.repute) == ("Dry earwax", "")
