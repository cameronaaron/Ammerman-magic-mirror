from __future__ import annotations

import pytest

from magic_mirror.genetics import (
    coerce_magnitude,
    normalize_genotype,
    parse_genotype_title,
)


def test_normalize_genotype_is_order_independent():
    assert normalize_genotype("AG") == normalize_genotype("GA") == "AG"
    assert normalize_genotype("gg") == "GG"


@pytest.mark.parametrize("bad", ["--", "DI", "A", "", "AN", "AGG"])
def test_normalize_genotype_rejects_non_diploid(bad):
    assert normalize_genotype(bad) is None


def test_parse_genotype_title_rsid_and_internal_id():
    assert parse_genotype_title("Rs17822931(T;T)") == ("rs17822931", "TT")
    assert parse_genotype_title("I3002432(A;A)") == ("i3002432", "AA")
    assert parse_genotype_title("Rs1(G;A)") == ("rs1", "AG")  # order-independent


@pytest.mark.parametrize("bad", ["Rs1(D;I)", "Rs1(A)", "Foo(A;A)", "", "rs1"])
def test_parse_genotype_title_rejects_malformed(bad):
    assert parse_genotype_title(bad) is None


@pytest.mark.parametrize(
    "value,expected",
    [("2.5", 2.5), (3, 3.0), ([4], 4.0), (["2.1"], 2.1), (None, 0.0), ("x", 0.0), ([], 0.0)],
)
def test_coerce_magnitude(value, expected):
    assert coerce_magnitude(value) == expected
