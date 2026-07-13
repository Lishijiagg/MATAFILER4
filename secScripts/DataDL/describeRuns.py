#!/usr/bin/env python3
"""
describeRuns.py - fetch descriptive ENA metadata for a set of run accessions.

Used to answer "what ARE these runs?" after an audit flags runs that a mapping
file does not use, or runs that expose no FASTQ.

Usage:
    ./describeRuns.py PRJNA703757 --missing-fastq
    ./describeRuns.py PRJNA675301 --runs SRR13004630 SRR13004631 ...
    ./describeRuns.py PRJEB45944  --not-in-map /path/to/PRJEB45944.map
"""

import argparse
import csv
import io
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ("run_accession,sample_accession,sample_title,sample_alias,"
          "library_name,library_strategy,library_layout,fastq_ftp,"
          "submitted_ftp,read_count,base_count,scientific_name,"
          "first_public,last_updated,study_title")
ACC_RE = re.compile(r"([EDS]RR\d+)")


def fetch(accession):
    params = urllib.parse.urlencode(
        {"accession": accession, "result": "read_run",
         "fields": FIELDS, "format": "tsv"})
    with urllib.request.urlopen(f"{ENA_PORTAL}?{params}", timeout=90) as r:
        tsv = r.read().decode("utf-8", "replace")
    return list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))


def map_accs(path):
    out = set()
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        found = ACC_RE.findall(cols[1] if len(cols) > 1 else cols[0])
        if found:
            out.add(max(found, key=len))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project")
    ap.add_argument("--runs", nargs="*", default=None,
                    help="only describe these run accessions")
    ap.add_argument("--not-in-map", type=Path, default=None,
                    help="describe the runs absent from this mapping file")
    ap.add_argument("--missing-fastq", action="store_true",
                    help="describe the runs that expose no fastq_ftp")
    args = ap.parse_args()

    rows = fetch(args.project)
    print(f"{args.project}: {len(rows)} run(s) in ENA\n")

    keep = rows
    if args.runs:
        want = set(args.runs)
        keep = [r for r in rows if r["run_accession"] in want]
    elif args.not_in_map:
        used = map_accs(args.not_in_map)
        keep = [r for r in rows if r["run_accession"] not in used]
        print(f"{len(keep)} run(s) absent from {args.not_in_map.name}\n")
    elif args.missing_fastq:
        keep = [r for r in rows if not (r.get("fastq_ftp") or "").strip()]
        print(f"{len(keep)} run(s) with no fastq_ftp\n")

    if not keep:
        print("nothing matched.")
        return 0

    for r in keep:
        has_fq = "yes" if (r.get("fastq_ftp") or "").strip() else "NO"
        has_sub = "yes" if (r.get("submitted_ftp") or "").strip() else "no"
        print(f"{r['run_accession']}")
        print(f"    sample       {r.get('sample_accession','')}")
        print(f"    title        {r.get('sample_title','')}")
        print(f"    alias        {r.get('sample_alias','')}")
        print(f"    library      {r.get('library_name','')}  "
              f"[{r.get('library_strategy','')}, {r.get('library_layout','')}]")
        print(f"    organism     {r.get('scientific_name','')}")
        print(f"    reads        {r.get('read_count','')}  "
              f"bases {r.get('base_count','')}")
        print(f"    fastq_ftp    {has_fq}      submitted_ftp {has_sub}")
        print(f"    first_public {r.get('first_public','')}   "
              f"last_updated {r.get('last_updated','')}")
        print()

    # a first_public date later than the rest of the study usually means the
    # runs were added after the data was originally downloaded
    dates = sorted({r.get("first_public", "") for r in rows if r.get("first_public")})
    kept_dates = sorted({r.get("first_public", "") for r in keep if r.get("first_public")})
    if dates and kept_dates and kept_dates[0] > dates[0]:
        print(f"NOTE: the whole study starts at first_public {dates[0]}, but "
              f"these runs start at {kept_dates[0]}. They may have been "
              f"deposited AFTER you downloaded the project.")
    return 0


if __name__ == "__main__":
    sys.exit(main())