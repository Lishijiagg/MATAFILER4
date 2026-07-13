#!/usr/bin/env python3
"""
auditMapIdentity.py - compare the run accessions in each mapping file against
what ENA actually holds, per library_strategy.

Counting runs is not enough: "58 in the map, 58 shotgun in ENA" does not prove
they are the *same* 58. This script compares the accession sets.

For every <ACCESSION>.map it reports:
    in map, not in ENA   -> the map references runs ENA does not return.
                            Never benign. Typos, duplicates, withdrawn runs,
                            or samples borrowed from another study.
    in ENA, not in map   -> runs you did not use, broken down by
                            library_strategy. Amplicon here is expected and
                            fine. Shotgun here needs an explanation.

Usage:
    ./auditMapIdentity.py <mapfile_dir> [--skip DIME] [--out report.tsv]
"""

import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,library_strategy,library_source"
ACC_RE = re.compile(r"([EDS]RR\d+)")


def ena_runs(accession, retries=4):
    """Return {run_accession: library_strategy} for a study."""
    params = urllib.parse.urlencode(
        {"accession": accession, "result": "read_run",
         "fields": FIELDS, "format": "tsv"})
    url = f"{ENA_PORTAL}?{params}"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                tsv = r.read().decode("utf-8", "replace")
            break
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            if attempt == retries:
                raise RuntimeError(f"ENA unreachable for {accession}: {exc}")
            time.sleep(min(30, 2 ** attempt))

    out = {}
    for row in csv.DictReader(io.StringIO(tsv), delimiter="\t"):
        acc = (row.get("run_accession") or "").strip()
        if acc:
            out[acc] = (row.get("library_strategy") or "NA").strip()
    return out


def map_accessions(path):
    """Pull run accessions out of a MATAFILER mapping file.

    The SmplPrefix column is the filename stem, e.g. 'ERR2587764_'. Take the
    longest accession-looking token in it. Returns {accession: [SmplID, ...]}
    so duplicate accessions become visible.
    """
    hits = defaultdict(list)
    unparsed = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            smpl_id = cols[0]
            prefix = cols[1] if len(cols) > 1 else ""
            found = ACC_RE.findall(prefix) or ACC_RE.findall(smpl_id)
            if found:
                hits[max(found, key=len)].append(smpl_id)
            else:
                unparsed.append(smpl_id)
    return hits, unparsed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mapdir", type=Path)
    ap.add_argument("--skip", nargs="*", default=["DIME"],
                    help="map stems that are not public accessions")
    ap.add_argument("--out", type=Path, default=None,
                    help="write a per-run TSV report here")
    args = ap.parse_args()

    report = []
    problems = []

    for m in sorted(args.mapdir.glob("*.map")):
        acc = m.stem
        if acc in args.skip:
            continue

        print(f"\n{'=' * 66}\n{acc}\n{'=' * 66}")
        in_map, unparsed = map_accessions(m)

        try:
            in_ena = ena_runs(acc)
        except RuntimeError as exc:
            print(f"  ENA ERROR: {exc}")
            problems.append((acc, "ena_unreachable", ""))
            continue

        dupes = {a: ids for a, ids in in_map.items() if len(ids) > 1}
        map_set, ena_set = set(in_map), set(in_ena)

        only_map = sorted(map_set - ena_set)
        only_ena = sorted(ena_set - map_set)

        print(f"  map: {sum(len(v) for v in in_map.values())} row(s), "
              f"{len(map_set)} distinct accession(s)")
        print(f"  ENA: {len(ena_set)} run(s)")

        if unparsed:
            print(f"  !! {len(unparsed)} map row(s) with no parsable "
                  f"accession: {', '.join(unparsed[:5])}")
            problems.append((acc, "unparsable_rows", ";".join(unparsed)))

        if dupes:
            print(f"  !! {len(dupes)} accession(s) appear on more than one "
                  f"map row:")
            for a, ids in list(dupes.items())[:5]:
                print(f"       {a} -> {', '.join(ids)}")
            problems.append((acc, "duplicate_accession",
                             ";".join(dupes)))

        if only_map:
            print(f"  !! {len(only_map)} accession(s) in the map but NOT in "
                  f"ENA -- this is never benign:")
            for a in only_map[:10]:
                print(f"       {a}  (SmplID: {', '.join(in_map[a])})")
            problems.append((acc, "in_map_not_in_ena", ";".join(only_map)))

        if only_ena:
            by_strat = defaultdict(list)
            for a in only_ena:
                by_strat[in_ena[a]].append(a)
            print(f"  {len(only_ena)} run(s) in ENA but not used:")
            for strat, accs in sorted(by_strat.items()):
                benign = strat.upper() in ("AMPLICON", "RNA-SEQ", "MIRNA-SEQ")
                tag = "expected" if benign else "** NEEDS EXPLANATION **"
                print(f"       {strat:<12} {len(accs):>4}   {tag}")
                if not benign:
                    print(f"         {', '.join(accs[:10])}"
                          f"{' ...' if len(accs) > 10 else ''}")
                    problems.append((acc, f"unused_{strat}", ";".join(accs)))

        if not (only_map or only_ena or dupes or unparsed):
            print("  OK: the map uses exactly the runs ENA returns.")

        for a in sorted(map_set | ena_set):
            report.append({
                "project": acc,
                "run": a,
                "in_map": "yes" if a in map_set else "no",
                "in_ena": "yes" if a in ena_set else "no",
                "library_strategy": in_ena.get(a, "NA"),
                "smpl_ids": ";".join(in_map.get(a, [])),
            })

    print(f"\n{'=' * 66}\nSUMMARY\n{'=' * 66}")
    if problems:
        for proj, kind, _ in problems:
            print(f"  {proj:<15} {kind}")
    else:
        print("  No discrepancies found.")

    if args.out and report:
        with args.out.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(report[0]), delimiter="\t")
            w.writeheader()
            w.writerows(report)
        print(f"\nPer-run report: {args.out}")

    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())