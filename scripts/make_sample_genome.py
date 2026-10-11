# SPDX-FileCopyrightText: 2026 Cameron Aaron
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Write a synthetic 23andMe-format genome for demos and tests.

Marker ids, chromosomes and positions come from a reference export (they are
the same for every customer of a given chip); every genotype is drawn at
random, so the output is nobody's DNA. A share of the markers SNPedia
describes get one of its low-magnitude genotypes, so the demo shows everyday
traits rather than a list of serious conditions.

    python scripts/make_sample_genome.py REFERENCE.txt result.csv sample/genome.txt
"""

import csv
import random
import re
import sys
from collections import defaultdict

EXTRA_MARKERS = 40_000
# Only these get a documented genotype, and only this often; the rest of the
# SNPedia markers get a random call like any other marker.
MAX_MAGNITUDE = 4.0
DOCUMENTED_SHARE = 0.25
CALLS = ["AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT"]


def snpedia_genotypes(path):
    """({id: genotypes safe to show}, {id: genotypes never to draw}).

    Never drawn: anything SNPedia marks Bad, or at magnitude MAX_MAGNITUDE
    and above. Ids are rsNNN or 23andMe's internal iNNN.
    """
    shown, never = defaultdict(list), defaultdict(set)
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.reader(f):
            title = row[0] if row else ""
            m = re.fullmatch(r"(Rs|I|rs|i)(\d+)\(([ACGT]);([ACGT])\)", title)
            if not m or len(row) < 3:
                continue
            marker = m.group(1).lower() + m.group(2)
            call = "".join(sorted(m.group(3) + m.group(4)))
            if row[2].strip() == "Bad" or float(row[1] or 0) >= MAX_MAGNITUDE:
                never[marker].add(call)
            else:
                shown[marker].append(call)
    return shown, never


def main(reference, snapshot, out):
    rng = random.Random(23)
    known, serious = snpedia_genotypes(snapshot)
    markers = []
    with open(reference, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rsid, chrom, pos = line.split("\t")[:3]
            markers.append((rsid, chrom, pos))
    documented_ids = set(known) | set(serious)
    chosen = [m for m in markers if m[0] in documented_ids]
    others = [m for m in markers if m[0] not in documented_ids]
    chosen += rng.sample(others, min(EXTRA_MARKERS, len(others)))
    order = {m: i for i, m in enumerate(markers)}
    chosen.sort(key=order.__getitem__)
    with open(out, "w", encoding="utf-8") as f:
        f.write("# Synthetic genome written by scripts/make_sample_genome.py.\n"
                "# Every genotype is random; this is not a real person's DNA.\n"
                "# rsid\tchromosome\tposition\tgenotype\n")
        for rsid, chrom, pos in chosen:
            if rsid in known and rng.random() < DOCUMENTED_SHARE:
                call = rng.choice(known[rsid])
            else:
                # Uniform calls are far more varied than a real genome and
                # would hit rare high-magnitude genotypes; never draw those.
                call = rng.choice([c for c in CALLS if c not in serious.get(rsid, ())]
                                  or CALLS)
            f.write(f"{rsid}\t{chrom}\t{pos}\t{call}\n")
    print(f"wrote {len(chosen)} markers to {out}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
