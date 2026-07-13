#!/usr/bin/env python3
"""
map2sampleTable.py - derive a fetchPublicReads sample table from an existing
MATAFILER mapping file.

A mapping file already encodes the study design (SmplID and AssmblGrps) that
fetchPublicReads cannot infer from ENA. If you have a curated map, this turns
it back into the table fetchPublicReads consumes, so a project can be
re-fetched reproducibly without redoing the curation by hand.

Usage:
    ./map2sampleTable.py PRJNA703757.map -o PRJNA703757.sample_table.tsv
"""

import argparse
import re
import sys
from pathlib import Path

ACC_RE = re.compile(r"([EDS]RR\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mapfile", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True)
    args = ap.parse_args()

    rows, skipped = [], []
    for line in args.mapfile.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        smpl_id = cols[0].strip()
        prefix = cols[1].strip() if len(cols) > 1 else ""
        grp = cols[2].strip() if len(cols) > 2 else ""

        found = ACC_RE.findall(prefix) or ACC_RE.findall(smpl_id)
        if not found:
            skipped.append(smpl_id)
            continue
        rows.append((max(found, key=len), smpl_id, grp))

    with args.out.open("w") as fh:
        fh.write("run_accession\tSmplID\tAssmblGrps\n")
        for acc, smpl_id, grp in rows:
            fh.write(f"{acc}\t{smpl_id}\t{grp}\n")

    print(f"{args.mapfile.name}: {len(rows)} sample(s) -> {args.out}")
    if skipped:
        print(f"  {len(skipped)} row(s) had no parsable run accession and "
              f"were skipped: {', '.join(skipped[:5])}", file=sys.stderr)

    groups = {g for _, _, g in rows if g}
    print(f"  {len(groups)} assembly group(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())