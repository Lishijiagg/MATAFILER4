#!/usr/bin/env python3
"""
fetchPublicReads.py - retrieve public metagenomic reads from ENA/SRA.

Part of MATAFILER4. Resolves run accessions via the ENA portal API,
downloads FASTQ files with wget (resumable), verifies them by MD5 and
writes a MATAFILER-compatible mapping file for the downloaded runs.

This script is a plain executable: it does NOT submit itself to a
scheduler. Job submission is handled by the calling MATAFILER routine
(Mods::GenoMetaAss::qsubSystem), so that cluster configuration stays in
one place.

Typical use:
    fetchPublicReads.py -a PRJNA755720 -o /scratch/raw/PRJNA755720 \
        --map-out /projects/meta/maps/PRJNA755720.map

    fetchPublicReads.py -m PRJNA755720_metadata.tsv -o /scratch/raw/PRJNA755720

Exit codes:
    0  all expected files present and MD5-verified
    1  incomplete: missing files, failed downloads or MD5 mismatches
    2  usage / input error
"""

import argparse
import csv
import hashlib
import io
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Fields requested from ENA. run_accession/fastq_ftp/fastq_md5 are required;
# the remainder are carried through into the mapping file / run report.
ENA_FIELDS = (
    "run_accession,sample_accession,secondary_sample_accession,"
    "study_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout,"
    "library_strategy,library_source,submitted_ftp,submitted_md5,"
    "instrument_platform,instrument_model,read_count,scientific_name"
)

# Column aliases accepted in user-supplied metadata tables. Keys are
# normalised (lower-case, stripped of spaces/underscores/hyphens).
URL_ALIASES = ("fastqftp", "fastqurl", "fastqaspera", "submittedftp")
MD5_ALIASES = ("fastqmd5", "submittedmd5")
RUN_ALIASES = ("runaccession", "accession", "run")
LAYOUT_ALIASES = ("librarylayout", "layout")
STRATEGY_ALIASES = ("librarystrategy", "strategy")
SOURCE_ALIASES = ("librarysource", "source")


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def log(msg: str) -> None:
    print(msg, flush=True)


def warn(msg: str) -> None:
    print("WARNING: " + msg, file=sys.stderr, flush=True)


def die(msg: str, code: int = 2) -> None:
    print("ERROR: " + msg, file=sys.stderr, flush=True)
    sys.exit(code)


def norm_key(s: str) -> str:
    return "".join(ch for ch in s.lower() if ch not in " _-")


def to_https(url: str) -> str:
    """ENA advertises bare ftp.sra.ebi.ac.uk paths; HTTPS is markedly more
    reliable for large transfers on HPC nodes behind a proxy."""
    url = url.strip()
    if not url:
        return url
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if not url.startswith(("http://", "https://")):
        return "https://" + url
    return url


# --------------------------------------------------------------------------
# data model
# --------------------------------------------------------------------------

@dataclass
class RemoteFile:
    url: str
    md5: Optional[str] = None
    size: int = -1        # ENA fastq_bytes; -1 = unknown

    @property
    def name(self) -> str:
        return os.path.basename(urllib.parse.urlparse(self.url).path)


@dataclass
class Run:
    accession: str
    files: List[RemoteFile] = field(default_factory=list)
    layout: str = ""
    platform: str = ""
    sample: str = ""
    study: str = ""
    strategy: str = ""   # ENA library_strategy: WGS, AMPLICON, RNA-Seq, ...
    source: str = ""     # ENA library_source: METAGENOMIC, GENOMIC, ...
    submitted: bool = False  # True if we fell back to submitted_ftp
    reads: int = -1      # ENA read_count; -1 = not reported

    @property
    def is_paired(self) -> bool:
        if self.layout:
            return self.layout.upper() == "PAIRED"
        return len([f for f in self.files if "_1." in f.name or "_2." in f.name]) >= 2


# --------------------------------------------------------------------------
# ENA portal API
# --------------------------------------------------------------------------

def ena_filereport(accession: str, retries: int = 4, timeout: int = 60) -> str:
    """Fetch the ENA filereport TSV for a study/sample/run accession."""
    params = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "tsv",
        }
    )
    url = f"{ENA_PORTAL}?{params}"
    last = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            if attempt < retries:
                wait = min(60, 2 ** attempt)
                warn(f"ENA query failed for {accession} "
                     f"(attempt {attempt}/{retries}): {exc}; retrying in {wait}s")
                time.sleep(wait)
    warn(f"ENA query permanently failed for {accession}: {last}")
    return ""


def parse_filereport(tsv: str, quiet: bool = False) -> List[Run]:
    """Turn an ENA filereport TSV into Run objects.

    Runs that carry no fastq_ftp are NOT silently dropped: ENA sometimes holds
    a run whose FASTQ files it never generated (only submitted_ftp exists, or
    processing never completed). Silently skipping those means quietly losing a
    sample, which is exactly the kind of error nobody notices. We warn, and
    fall back to submitted_ftp when it is available.
    """
    runs: List[Run] = []
    no_files: List[str] = []
    fallback: List[str] = []
    if not tsv.strip():
        return runs
    reader = csv.DictReader(io.StringIO(tsv), delimiter="\t")
    for row in reader:
        acc = (row.get("run_accession") or "").strip()
        if not acc:
            continue

        ftp = (row.get("fastq_ftp") or "").strip()
        md5_field = (row.get("fastq_md5") or "").strip()
        used_fallback = False

        if not ftp:
            # fall back to the submitter's own upload
            ftp = (row.get("submitted_ftp") or "").strip()
            md5_field = (row.get("submitted_md5") or "").strip()
            used_fallback = bool(ftp)
            if not ftp:
                no_files.append(acc)
                continue
            fallback.append(acc)

        urls = [u for u in ftp.split(";") if u.strip()]
        md5s = [m.strip() for m in md5_field.split(";")]
        if len(md5s) != len(urls):
            md5s = [None] * len(urls)
        raw_sizes = [b.strip() for b in (row.get("fastq_bytes") or "").split(";")]
        sizes = ([int(b) if b.isdigit() else -1 for b in raw_sizes]
                 if len(raw_sizes) == len(urls) else [-1] * len(urls))
        run = Run(
            accession=acc,
            layout=(row.get("library_layout") or "").strip(),
            strategy=(row.get("library_strategy") or "").strip(),
            source=(row.get("library_source") or "").strip(),
            platform=(row.get("instrument_platform") or "").strip(),
            sample=(row.get("sample_accession") or "").strip(),
            study=(row.get("study_accession") or "").strip(),
            submitted=used_fallback,
            reads=int(row["read_count"]) if (row.get("read_count") or "").strip().isdigit() else -1,
            files=[RemoteFile(to_https(u), m or None, sz)
                   for u, m, sz in zip(urls, md5s, sizes)],
        )
        runs.append(run)

    if not quiet and fallback:
        warn(f"{len(fallback)} run(s) have no ENA-generated FASTQ; falling "
             f"back to submitted_ftp: {', '.join(fallback[:5])}"
             f"{' ...' if len(fallback) > 5 else ''}. These files are as the "
             f"submitter uploaded them -- check the format before use.")
    if not quiet and no_files:
        warn(f"{len(no_files)} run(s) exist in ENA but expose NO downloadable "
             f"files at all (neither fastq_ftp nor submitted_ftp) and were "
             f"skipped: {', '.join(no_files[:5])}"
             f"{' ...' if len(no_files) > 5 else ''}. If any of these belong "
             f"in your study you must obtain them another way (SRA).")
    return runs


# --------------------------------------------------------------------------
# user-supplied metadata table
# --------------------------------------------------------------------------

def parse_metadata(path: Path) -> List[Run]:
    """Parse a user-supplied metadata TSV.

    Accepts both ENA-native headers (fastq_ftp / fastq_md5 / run_accession)
    and the ENA browser export style ('FastQ URL', 'FastQ MD5', 'Accession').
    Header matching ignores case, spaces, underscores and hyphens. A UTF-8
    BOM, if present, is stripped.
    """
    runs: List[Run] = []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            die(f"metadata file has no header: {path}")

        lookup = {norm_key(f): f for f in reader.fieldnames if f}
        url_col = next((lookup[a] for a in URL_ALIASES if a in lookup), None)
        md5_col = next((lookup[a] for a in MD5_ALIASES if a in lookup), None)
        run_col = next((lookup[a] for a in RUN_ALIASES if a in lookup), None)
        lay_col = next((lookup[a] for a in LAYOUT_ALIASES if a in lookup), None)
        strat_col = next((lookup[a] for a in STRATEGY_ALIASES if a in lookup), None)
        src_col = next((lookup[a] for a in SOURCE_ALIASES if a in lookup), None)

        if not url_col:
            die(
                f"no FASTQ URL column found in {path}.\n"
                f"       looked for (case/space insensitive): {', '.join(URL_ALIASES)}\n"
                f"       header was: {', '.join(reader.fieldnames)}"
            )
        log(f"  URL column      : {url_col}")
        log(f"  MD5 column      : {md5_col or '<absent, will query ENA>'}")
        log(f"  accession column: {run_col or '<absent, will infer from URL>'}")

        for row in reader:
            raw = (row.get(url_col) or "").strip()
            if not raw:
                continue
            urls = [u.strip() for u in raw.split(";") if u.strip()]
            md5s = [m.strip() for m in (row.get(md5_col) or "").split(";")] if md5_col else []
            if len(md5s) != len(urls):
                md5s = [None] * len(urls)

            acc = (row.get(run_col) or "").strip() if run_col else ""
            if not acc:
                acc = infer_accession(urls[0])
            if not acc:
                warn(f"cannot determine run accession for {urls[0]}, skipping")
                continue

            runs.append(
                Run(
                    accession=acc,
                    layout=(row.get(lay_col) or "").strip() if lay_col else "",
                    strategy=(row.get(strat_col) or "").strip() if strat_col else "",
                    source=(row.get(src_col) or "").strip() if src_col else "",
                    files=[RemoteFile(to_https(u), m or None) for u, m in zip(urls, md5s)],
                )
            )
    return runs


ACC_RE = re.compile(r"([EDS]RR\d+)")


def infer_accession(url: str) -> str:
    """Pull an (E|S|D)RR accession out of a URL.

    Careful: ENA paths embed *truncated* accessions as directory names,
    e.g. /vol1/fastq/SRR154/051/SRR15423451/SRR15423451_1.fastq.gz
                     ^^^^^^ not the accession
    Match against the basename only, and take the longest hit.
    """
    base = os.path.basename(urllib.parse.urlparse(url).path)
    hits = ACC_RE.findall(base) or ACC_RE.findall(url)
    return max(hits, key=len) if hits else ""


# --------------------------------------------------------------------------
# MD5 backfill
# --------------------------------------------------------------------------

def backfill_md5(runs: Sequence[Run], threads: int) -> None:
    """Query ENA for MD5 sums of any run whose files lack them."""
    need = [r for r in runs if any(f.md5 is None for f in r.files)]
    if not need:
        return
    log(f"\n=== Backfilling MD5 from ENA for {len(need)} run(s) ===")

    def one(run: Run) -> None:
        fetched = parse_filereport(ena_filereport(run.accession), quiet=True)
        if not fetched:
            return
        table = {f.name: f.md5 for f in fetched[0].files if f.md5}
        for f in run.files:
            if f.md5 is None and f.name in table:
                f.md5 = table[f.name]

    with ThreadPoolExecutor(max_workers=min(threads, 8)) as pool:
        list(as_completed([pool.submit(one, r) for r in need]))

    still = sum(1 for r in runs for f in r.files if f.md5 is None)
    if still:
        warn(f"{still} file(s) still have no MD5 after ENA backfill")


# --------------------------------------------------------------------------
# download
# --------------------------------------------------------------------------

def wget_one(rf: RemoteFile, outdir: Path, wget_bin: str, opts: argparse.Namespace) -> bool:
    """Download a single file with wget. Returns True on success.

    wget (rather than a pure-Python client) is used deliberately: its
    resume (-c), retry and timeout handling is battle-tested against the
    ENA mirrors, which drop long transfers fairly regularly.
    """
    dest = outdir / rf.name

    # Deciding whether a file is already done is where this used to be wrong.
    # "exists and is non-empty" is NOT the same as "complete": an interrupted
    # download leaves a truncated file that is non-empty, and skipping it means
    # silently keeping corrupt data. ENA tells us the exact byte count, so use
    # it.
    if dest.is_file() and not opts.force:
        on_disk = dest.stat().st_size
        if rf.size > 0:
            if on_disk == rf.size:
                return True                       # genuinely complete
            if on_disk > rf.size:
                warn(f"{rf.name} is larger than ENA declares "
                     f"({on_disk} > {rf.size} bytes); discarding and "
                     f"re-fetching")
                dest.unlink()
            # on_disk < rf.size: leave it, wget -c resumes from there
        elif on_disk > 0:
            # ENA gave us no size. We cannot tell complete from truncated, so
            # let wget -c try to resume; MD5 verification is the safety net.
            pass

    cmd = [
        wget_bin,
        "-c",
        "--tries", str(opts.retries),
        "--waitretry", str(opts.waitretry),
        "--connect-timeout", str(opts.connect_timeout),
        "--read-timeout", str(opts.read_timeout),
        "--no-verbose",
        "-P", str(outdir),      # NOT -O: wget warns that -O breaks -c resume
        rf.url,
    ]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except OSError as exc:
        warn(f"cannot execute {wget_bin}: {exc}")
        return False
    if res.returncode != 0:
        warn(f"download failed ({rf.url}): "
             f"{res.stdout.decode('utf-8', 'replace').strip()[-200:]}")
        return False
    if not dest.is_file() or dest.stat().st_size == 0:
        return False
    if dest.name.endswith(".gz") and not looks_like_fastq_gz(dest):
        head = dest.read_bytes()[:80].decode("utf-8", "replace").replace("\n", " ")
        warn(f"{rf.name} is not gzip -- ENA likely returned an error page "
             f"(starts: {head!r}). Discarding.")
        dest.unlink()
        return False
    if rf.size > 0 and dest.stat().st_size != rf.size:
        warn(f"{rf.name} is {dest.stat().st_size} bytes, ENA declares "
             f"{rf.size}. Treating as failed.")
        return False
    return True


GZIP_MAGIC = b"\x1f\x8b"


def looks_like_fastq_gz(path: Path) -> bool:
    """Guard against ENA handing back an HTML error page (a directory-index
    listing, a 404, a rate-limit notice) with a .fastq.gz name. Real gzip
    starts with the two magic bytes 1f 8b; an HTML page starts with '<' or
    whitespace. Also reject anything implausibly small: a real metagenomic
    FASTQ is megabytes at minimum, so a sub-kilobyte .gz is always junk."""
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with path.open("rb") as fh:
            magic = fh.read(2)
    except OSError:
        return False
    # the definitive test is the gzip magic number; an HTML error page fails it
    return magic == GZIP_MAGIC


def complete_on_disk(rf: RemoteFile, outdir: Path) -> bool:
    """A file counts as complete only if its size matches ENA's fastq_bytes.

    Falling back to "non-empty" when ENA reports no size is the best we can do
    there; MD5 verification catches the rest.
    """
    dest = outdir / rf.name
    if not dest.is_file():
        return False
    # an HTML error page saved as .fastq.gz has the "right" name and a nonzero
    # size but is not gzip; reject it regardless of the byte count
    if dest.name.endswith((".gz",)) and not looks_like_fastq_gz(dest):
        return False
    on_disk = dest.stat().st_size
    if rf.size > 0:
        return on_disk == rf.size
    return on_disk > 0


def download_all(runs: Sequence[Run], outdir: Path, wget_bin: str,
                 opts: argparse.Namespace) -> List[RemoteFile]:
    """Two-pass download. Returns the files that are still missing."""
    files = [f for r in runs for f in r.files]
    log(f"\n=== Pass 1: downloading {len(files)} file(s), "
        f"{opts.threads} parallel stream(s) ===")

    started = time.time()
    with ThreadPoolExecutor(max_workers=opts.threads) as pool:
        futures = {pool.submit(wget_one, f, outdir, wget_bin, opts): f for f in files}
        done = 0
        for fut in as_completed(futures):
            rf = futures[fut]
            done += 1
            dest = outdir / rf.name
            size = dest.stat().st_size if dest.is_file() else 0
            mins = (time.time() - started) / 60.0
            log(f"  [{done}/{len(files)}] {rf.name}  "
                f"{size / 1e9:.2f} GB  ({mins:.1f} min elapsed)")

    missing = [f for f in files if not complete_on_disk(f, outdir)]
    log(f"Pass 1 complete. Missing or incomplete: {len(missing)}")
    if not missing:
        return []

    # Pass 2+: ENA failures tend to come in windows lasting minutes to tens of
    # minutes -- a mirror node hiccups and every retry inside that window hits
    # the same wall. wget's own --tries burns through in seconds, well inside
    # the window, so retrying immediately is nearly useless. We therefore wait
    # a growing amount of time between passes to climb out of the window, and
    # crucially re-check completeness with complete_on_disk (NOT size>0): a
    # truncated file is non-empty and must be retried, not skipped.
    for attempt in range(1, opts.resume_passes + 1):
        wait = min(opts.resume_wait * attempt, 900)
        log(f"\n=== Pass {attempt + 1}: {len(missing)} file(s) still "
            f"incomplete; waiting {wait}s for ENA to settle, then retrying ===")
        time.sleep(wait)

        stale_runs = sorted({r.accession for r in runs
                             for f in r.files if f in missing})
        for acc in stale_runs:
            fresh = parse_filereport(ena_filereport(acc), quiet=True)
            if not fresh:
                warn(f"ENA returned nothing for {acc}")
                continue
            for rf in fresh[0].files:
                # only touch files that are still not complete
                if complete_on_disk(rf, outdir):
                    continue
                log(f"  retrying {rf.name}")
                if wget_one(rf, outdir, wget_bin, opts):
                    for r in runs:
                        for f in r.files:
                            if f.name == rf.name:
                                if rf.md5:
                                    f.md5 = rf.md5
                                if rf.size > 0:
                                    f.size = rf.size

        missing = [f for f in files if not complete_on_disk(f, outdir)]
        log(f"Pass {attempt + 1} complete. Still incomplete: {len(missing)}")
        if not missing:
            return []

    return missing


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------

def md5_of(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def verify_all(runs: Sequence[Run], outdir: Path, threads: int) -> List[str]:
    """MD5-verify every downloaded file. Returns a list of failure messages.

    A matching MD5 also proves gzip-stream integrity, so no separate
    `gzip -t` pass is needed.

    This pass is pure I/O and can take a long time -- hashing several hundred
    gigabytes off a shared filesystem is easily half an hour -- so it reports
    as it goes. A silent half hour is indistinguishable from a hang.
    """
    targets = [f for r in runs for f in r.files if f.md5]
    skipped = sum(1 for r in runs for f in r.files if not f.md5)
    total_bytes = sum((outdir / f.name).stat().st_size
                      for f in targets if (outdir / f.name).is_file())
    log(f"\n=== Verifying MD5 for {len(targets)} file(s), "
        f"{total_bytes / 1e9:.0f} GB "
        f"({skipped} without reference MD5) ===")
    log(f"    This is I/O bound and reports as it goes. Be patient.")

    failures: List[str] = []
    started = time.time()

    def one(rf: RemoteFile) -> Optional[str]:
        path = outdir / rf.name
        if not path.is_file() or path.stat().st_size == 0:
            return f"MISSING  {rf.name}"
        actual = md5_of(path)
        if actual != rf.md5:
            return f"MISMATCH {rf.name} (expected {rf.md5}, got {actual})"
        return None

    done = 0
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(one, f): f for f in targets}
        for fut in as_completed(futures):
            done += 1
            msg = fut.result()
            if msg:
                failures.append(msg)
                log(f"  [{done}/{len(targets)}] {msg}")
            elif done % 25 == 0 or done == len(targets):
                mins = (time.time() - started) / 60.0
                rate = done / mins if mins > 0 else 0
                eta = (len(targets) - done) / rate if rate > 0 else 0
                log(f"  [{done}/{len(targets)}] ok  "
                    f"({mins:.1f} min elapsed, ~{eta:.0f} min left)")

    log(f"MD5 verification done in {(time.time() - started) / 60.0:.1f} min. "
        f"Failures: {len(failures)}")
    return failures


# --------------------------------------------------------------------------
# run selection
#
# A study accession frequently resolves to more runs than you want. PRJEB26842
# holds 116 runs: 58 amplicon and 58 shotgun runs of the same samples. Feeding
# amplicon reads into an assembly pipeline produces nonsense silently, so
# selection is explicit rather than best-effort.
#
# Note on library_strategy: it is declared by the submitter and is NOT
# reliable. The shotgun runs of PRJEB26842 are labelled OTHER, not WGS. So:
#
#   --sample-table          authoritative. Your table, your runs. Preferred.
#   --exclude-strategy      robust. AMPLICON is rarely mislabelled as anything
#                           else, so excluding it is safe across studies.
#   --library-strategy      fragile. Only use when you have checked what the
#                           submitter actually wrote.
#
# Order matters: the whitelist is applied first, so that later filters cannot
# produce misleading "not found in ENA" diagnostics.
# --------------------------------------------------------------------------

def select_runs(runs: List[Run], opts: argparse.Namespace,
                sample_table: Optional[Dict[str, Dict[str, str]]]) -> List[Run]:
    resolved = {r.accession for r in runs}
    kept = runs

    # 1. sample table whitelist -- authoritative, applied first
    if sample_table is not None and not opts.no_whitelist:
        before = len(kept)
        kept = [r for r in kept if r.accession in sample_table]
        log(f"  sample table whitelist: {before} -> {len(kept)} run(s)")
        absent = sorted(set(sample_table) - resolved)
        if absent:
            warn(f"{len(absent)} run(s) in the sample table are absent from "
                 f"ENA's response: {', '.join(absent[:5])}"
                 f"{' ...' if len(absent) > 5 else ''}")

    # 2. empty / near-empty runs
    #
    # ENA studies accumulate broken submissions. PRJNA703757 gained a run in
    # 2025 with read_count 0 and no files at all, four years after the study
    # was published. Such a run in a mapping file is a phantom sample: the
    # pipeline finds nothing for it and your sample count is quietly wrong.
    empty = [r for r in kept if r.reads == 0]
    if empty:
        warn(f"{len(empty)} run(s) report read_count 0 and were dropped: "
             f"{', '.join(r.accession for r in empty[:5])}"
             f"{' ...' if len(empty) > 5 else ''}")
        kept = [r for r in kept if r.reads != 0]

    if opts.min_reads > 0:
        before = len(kept)
        thin = [r for r in kept if 0 <= r.reads < opts.min_reads]
        kept = [r for r in kept if r.reads < 0 or r.reads >= opts.min_reads]
        log(f"  min-reads {opts.min_reads}: {before} -> {len(kept)} run(s)")
        for r in thin[:10]:
            log(f"    dropped {r.accession} ({r.reads} reads)")

    # 3. exclusion -- robust
    if opts.exclude_strategy:
        drop = {s.strip().upper() for s in opts.exclude_strategy.split(",")}
        before = len(kept)
        removed = [r for r in kept if r.strategy.upper() in drop]
        kept = [r for r in kept if r.strategy.upper() not in drop]
        log(f"  excluding library_strategy {sorted(drop)}: "
            f"{before} -> {len(kept)} run(s)")
        if removed:
            log(f"    dropped: {', '.join(r.accession for r in removed[:5])}"
                f"{' ...' if len(removed) > 5 else ''}")

    # 4. inclusion -- fragile, opt-in
    if opts.library_strategy:
        want = {s.strip().upper() for s in opts.library_strategy.split(",")}
        before = len(kept)
        blank = [r for r in kept if not r.strategy]
        kept = [r for r in kept if r.strategy.upper() in want]
        log(f"  keeping library_strategy {sorted(want)}: "
            f"{before} -> {len(kept)} run(s)")
        if blank:
            warn(f"{len(blank)} run(s) declare no library_strategy and were "
                 f"dropped")
        if before and not kept:
            warn("the library_strategy filter removed every run. "
                 "library_strategy is submitter-declared and often wrong "
                 "(PRJEB26842 labels its shotgun runs OTHER, not WGS). "
                 "Check the breakdown printed above, or prefer "
                 "--sample-table / --exclude-strategy.")

    if opts.library_source:
        want = {s.strip().upper() for s in opts.library_source.split(",")}
        before = len(kept)
        kept = [r for r in kept if r.source.upper() in want]
        log(f"  keeping library_source {sorted(want)}: "
            f"{before} -> {len(kept)} run(s)")

    return kept


def strategy_breakdown(runs: Sequence[Run]) -> str:
    counts: Dict[str, int] = {}
    for r in runs:
        key = f"{r.strategy or 'NA'}/{r.source or 'NA'}"
        counts[key] = counts.get(key, 0) + 1
    return ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))


def warn_if_heterogeneous(runs: Sequence[Run], opts: argparse.Namespace,
                          sample_table: Optional[Dict]) -> None:
    """A study carrying more than one library_strategy is the classic trap."""
    kinds = {r.strategy.upper() for r in runs if r.strategy}
    if len(kinds) < 2:
        return
    filtered = (opts.library_strategy or opts.exclude_strategy
                or (sample_table is not None and not opts.no_whitelist))
    if filtered:
        return
    warn(f"this study mixes {len(kinds)} library strategies "
         f"({', '.join(sorted(kinds))}) and NO run selection was requested. "
         f"All of them will be downloaded and written to the mapping file. "
         f"If some are amplicon runs, assembling them will silently produce "
         f"nonsense. Use --sample-table, --exclude-strategy or "
         f"--library-strategy.")


# --------------------------------------------------------------------------
# MATAFILER mapping file
#
# Two distinct classes of information end up in a mapping file:
#
#   Derivable here      SmplPrefix, #DirPath, #RunID
#                       -- we downloaded the files, so we know their names.
#
#   NOT derivable here  SmplID, AssmblGrps
#                       -- these encode experimental design (which runs are
#                          the same subject and should be co-assembled; what
#                          a sample is called biologically). ENA does not
#                          carry that reliably, and guessing it from sample
#                          titles is study-specific.
#
# The second class is therefore supplied by the user through --sample-table,
# a documented TSV keyed on run accession. Without it we fall back to safe,
# study-agnostic defaults: SmplID = run accession, and every sample in its
# own assembly group (i.e. no co-assembly).
# --------------------------------------------------------------------------

SAMPLE_TABLE_ID = ("smplid", "sampleid", "samplename", "sample")
SAMPLE_TABLE_GRP = ("assmblgrps", "assemblygroup", "assemblygroups", "group")
SAMPLE_TABLE_RUN = ("runaccession", "accession", "run")


def read_sample_table(path: Path) -> Dict[str, Dict[str, str]]:
    """Read the user's run -> (SmplID, AssmblGrps) table.

    Required column: a run accession. Optional: SmplID, AssmblGrps.
    Header matching ignores case, spaces, underscores and hyphens.
    """
    table: Dict[str, Dict[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            die(f"sample table has no header: {path}")
        lookup = {norm_key(f): f for f in reader.fieldnames if f}
        run_col = next((lookup[a] for a in SAMPLE_TABLE_RUN if a in lookup), None)
        id_col = next((lookup[a] for a in SAMPLE_TABLE_ID if a in lookup), None)
        grp_col = next((lookup[a] for a in SAMPLE_TABLE_GRP if a in lookup), None)
        if not run_col:
            die(f"no run accession column in {path}. "
                f"Expected one of: {', '.join(SAMPLE_TABLE_RUN)}")
        if not id_col and not grp_col:
            die(f"{path} carries neither a SmplID nor an AssmblGrps column; "
                f"it would add nothing")
        log(f"  sample table: run={run_col}, "
            f"SmplID={id_col or '<default>'}, AssmblGrps={grp_col or '<default>'}")
        for row in reader:
            acc = (row.get(run_col) or "").strip()
            if not acc:
                continue
            table[acc] = {
                "SmplID": (row.get(id_col) or "").strip() if id_col else "",
                "AssmblGrps": (row.get(grp_col) or "").strip() if grp_col else "",
            }
    return table


def smpl_prefix(run: Run) -> str:
    """Derive the SmplPrefix MATAFILER globs DirPath with.

    Taken from the filenames ENA declares, NOT from what is on disk, so the
    mapping file can be produced -- and reviewed -- before committing to a
    download that may run to hundreds of gigabytes.

    Paired runs give ERR123_1.fastq.gz + ERR123_2.fastq.gz -> "ERR123_".
    Single-end gives ERR123.fastq.gz -> "ERR123." (keeping the dot so the
    prefix cannot also match ERR1234.fastq.gz).
    """
    names = sorted(f.name for f in run.files)
    if not names:
        return ""
    if len(names) == 1:
        stem = names[0]
        for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if stem.endswith(suf):
                return stem[: -len(suf)] + "."
        return stem
    prefix = os.path.commonprefix(names)
    return prefix or (run.accession + "_")


def write_mapping(runs: Sequence[Run], outdir: Path, map_path: Path,
                  opts: argparse.Namespace,
                  sample_table: Optional[Dict[str, Dict[str, str]]],
                  check_disk: bool = True) -> int:
    if check_disk:
        ok = [r for r in runs
              if r.files and all((outdir / f.name).is_file()
                                 and (outdir / f.name).stat().st_size > 0
                                 for f in r.files)]
    else:
        ok = [r for r in runs if r.files]

    seen_ids: Dict[str, str] = {}
    rows: List[tuple] = []
    unannotated: List[str] = []

    for r in ok:
        entry = (sample_table or {}).get(r.accession, {})
        smpl_id = entry.get("SmplID") or r.accession
        grp = entry.get("AssmblGrps") or smpl_id  # own group == no co-assembly
        if sample_table is not None and r.accession not in sample_table:
            unannotated.append(r.accession)

        if smpl_id in seen_ids and seen_ids[smpl_id] != r.accession:
            die(f"duplicate SmplID '{smpl_id}' for runs {seen_ids[smpl_id]} "
                f"and {r.accession}. MATAFILER requires unique sample IDs; "
                f"fix the sample table.", 1)
        seen_ids[smpl_id] = r.accession

        prefix = smpl_prefix(r)
        if not prefix:
            warn(f"cannot derive SmplPrefix for {r.accession}, skipping")
            continue
        rows.append((smpl_id, prefix, grp))

    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("w") as fh:
        fh.write("#SmplID\tSmplPrefix\tAssmblGrps\n")
        if opts.out_path:
            fh.write(f"#OutPath\t{opts.out_path}\n")
        fh.write(f"#RunID\t{opts.run_id}\n")
        if opts.relax_smplid:
            fh.write("#RelaxSMPLID\tTRUE\n")
        fh.write(f"#DirPath\t{outdir}\n")
        for smpl_id, prefix, grp in rows:
            fh.write(f"{smpl_id}\t{prefix}\t{grp}\n")

    log(f"\nMapping file written: {map_path} ({len(rows)} sample(s))")
    if check_disk and len(ok) != len(runs):
        warn(f"{len(runs) - len(ok)} run(s) omitted: not all their files "
             f"are present on disk")

    # a co-assembly group holding a single sample is worth surfacing: it may
    # be intentional (a subject with one timepoint) or a sign that the rest of
    # the group was lost somewhere upstream
    counts: Dict[str, int] = {}
    for _, _, g in rows:
        counts[g] = counts.get(g, 0) + 1
    singles = [g for g, n in counts.items() if n == 1]
    if singles:
        log(f"  {len(counts)} assembly group(s); {len(singles)} contain a "
            f"single sample (no co-assembly for those)")
    if unannotated:
        warn(f"{len(unannotated)} run(s) absent from the sample table, "
             f"defaulted to accession / no co-assembly: "
             f"{', '.join(unannotated[:5])}"
             f"{' ...' if len(unannotated) > 5 else ''}")
    if sample_table is None:
        log("  NOTE: no --sample-table given. SmplID = run accession and every "
            "sample is its own assembly group (no co-assembly). Supply "
            "--sample-table to encode the study design.")
    return len(rows)


def write_report(runs: Sequence[Run], outdir: Path, report: Path) -> None:
    with report.open("w") as fh:
        fh.write("run_accession\tfile\tsize_bytes\tmd5_expected\tstatus\n")
        for r in runs:
            for f in r.files:
                p = outdir / f.name
                if p.is_file() and p.stat().st_size > 0:
                    status, size = "present", p.stat().st_size
                else:
                    status, size = "missing", 0
                fh.write(f"{r.accession}\t{f.name}\t{size}\t{f.md5 or 'NA'}\t{status}\n")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fetchPublicReads.py",
        description="Download public metagenomic reads from ENA/SRA and emit "
                    "a MATAFILER mapping file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession",
                     help="study/project accession to resolve via ENA "
                          "(e.g. PRJNA755720, PRJEB1234, ERP001234)")
    src.add_argument("-m", "--metadata", type=Path,
                     help="pre-downloaded metadata TSV (ENA filereport or "
                          "ENA browser export)")

    p.add_argument("-o", "--outdir", type=Path, required=True,
                   help="directory the FASTQ files are written to")
    p.add_argument("--map-out", type=Path, default=None,
                   help="path of the MATAFILER mapping file to write "
                        "[default: <outdir>/<accession>.map]")
    p.add_argument("--out-path", default=None,
                   help="value written to the #OutPath header of the mapping "
                        "file (MATAFILER's output directory)")
    p.add_argument("--run-id", default=None,
                   help="value written to the #RunID header "
                        "[default: the accession / metadata file stem]")
    p.add_argument("--sample-table", type=Path, default=None,
                   help="TSV supplying the study design, keyed on run "
                        "accession. Recognised columns: run_accession "
                        "(required), SmplID, AssmblGrps. Runs sharing an "
                        "AssmblGrps value are co-assembled. Without this, "
                        "SmplID defaults to the run accession and no samples "
                        "are co-assembled")
    p.add_argument("--relax-smplid", action="store_true",
                   help="write #RelaxSMPLID TRUE into the mapping file")

    p.add_argument("--exclude-strategy", default=None,
                   help="drop runs with these ENA library_strategy values "
                        "(comma-separated, e.g. AMPLICON). Robust: AMPLICON "
                        "is rarely mislabelled, so excluding it works across "
                        "studies. Prefer this over --library-strategy")
    p.add_argument("--library-strategy", default=None,
                   help="only keep runs with these ENA library_strategy "
                        "values (comma-separated). FRAGILE: library_strategy "
                        "is submitter-declared and often wrong -- PRJEB26842 "
                        "labels its shotgun runs OTHER, not WGS. Check the "
                        "breakdown the tool prints before relying on this")
    p.add_argument("--library-source", default=None,
                   help="only keep runs with these ENA library_source values "
                        "(comma-separated, e.g. METAGENOMIC)")
    p.add_argument("--min-reads", type=int, default=0,
                   help="drop runs with fewer than this many reads. Useful "
                        "for excluding blanks and failed libraries. Runs "
                        "reporting read_count 0 are ALWAYS dropped")
    p.add_argument("--no-whitelist", action="store_true",
                   help="do not treat --sample-table as a whitelist; keep "
                        "every run ENA returns and merely annotate those "
                        "present in the table")

    p.add_argument("-t", "--threads", type=int, default=4,
                   help="parallel download / verification streams. Keep this "
                        "low: ENA throttles aggressive clients")
    p.add_argument("--retries", type=int, default=20, help="wget --tries")
    p.add_argument("--waitretry", type=int, default=60, help="wget --waitretry")
    p.add_argument("--resume-passes", type=int, default=4,
                   help="how many extra retry passes after the first, each "
                        "preceded by a growing wait, to climb out of an ENA "
                        "outage window")
    p.add_argument("--resume-wait", type=int, default=120,
                   help="base seconds to wait before retry pass N; the actual "
                        "wait is resume-wait * N, capped at 900")
    p.add_argument("--connect-timeout", type=int, default=60)
    p.add_argument("--read-timeout", type=int, default=600)

    p.add_argument("--wget", dest="wget_bin", default=None,
                   help="path to wget [default: from PATH]")
    p.add_argument("--force", action="store_true",
                   help="re-download files that already exist")
    p.add_argument("--no-verify", action="store_true",
                   help="skip MD5 verification (not recommended)")
    p.add_argument("--dry-run", action="store_true",
                   help="resolve accessions and report what would be "
                        "downloaded, then exit")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    opts = build_parser().parse_args(argv)

    wget_bin = opts.wget_bin or shutil.which("wget")
    if not wget_bin:
        die("wget not found. Pass --wget or add it to PATH.")

    label = opts.accession or opts.metadata.stem
    if not opts.run_id:
        opts.run_id = label

    sample_table: Optional[Dict[str, Dict[str, str]]] = None
    if opts.sample_table:
        if not opts.sample_table.is_file():
            die(f"sample table not found: {opts.sample_table}")
        log(f"=== Reading sample table {opts.sample_table} ===")
        sample_table = read_sample_table(opts.sample_table)
        log(f"  {len(sample_table)} run(s) annotated")

    log("=" * 60)
    log(f"MATAFILER4 :: fetchPublicReads")
    log(f"Source     : {opts.accession or opts.metadata}")
    log(f"Output dir : {opts.outdir}")
    log(f"Started    : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 60)

    # 1. resolve
    if opts.accession:
        log(f"\n=== Resolving {opts.accession} via ENA portal ===")
        runs = parse_filereport(ena_filereport(opts.accession))
    else:
        if not opts.metadata.is_file() or opts.metadata.stat().st_size == 0:
            die(f"metadata file missing or empty: {opts.metadata}")
        log(f"\n=== Parsing {opts.metadata} ===")
        runs = parse_metadata(opts.metadata)

    if not runs:
        die(f"no runs with FASTQ files resolved for {label}", 1)

    log(f"Resolved {len(runs)} run(s) / "
        f"{sum(len(r.files) for r in runs)} file(s)")
    log(f"  library_strategy/source: {strategy_breakdown(runs)}")

    warn_if_heterogeneous(runs, opts, sample_table)

    log("\n=== Selecting runs ===")
    runs = select_runs(runs, opts, sample_table)
    if not runs:
        die("no runs left after filtering; check --library-strategy / "
            "--sample-table", 1)

    n_files = sum(len(r.files) for r in runs)
    log(f"\nSelected {len(runs)} run(s) / {n_files} file(s)")

    if opts.dry_run:
        map_path = opts.map_out or (opts.outdir / f"{label}.map")
        write_mapping(runs, opts.outdir, map_path, opts, sample_table,
                      check_disk=False)
        log("\nDRY RUN: nothing was downloaded. The mapping file above is "
            "what a real run would produce -- diff it against your existing "
            "map before committing to the download.")
        return 0

    opts.outdir.mkdir(parents=True, exist_ok=True)

    # 2. MD5 backfill
    if not opts.no_verify:
        backfill_md5(runs, opts.threads)

    # 3. download
    missing = download_all(runs, opts.outdir, wget_bin, opts)

    # 4. verify
    failures: List[str] = []
    if opts.no_verify:
        warn("MD5 verification skipped by request")
    else:
        failures = verify_all(runs, opts.outdir, opts.threads)

    # 5. outputs
    map_path = opts.map_out or (opts.outdir / f"{label}.map")
    write_mapping(runs, opts.outdir, map_path, opts, sample_table)
    report = opts.outdir / f"{label}.download_report.tsv"
    write_report(runs, opts.outdir, report)

    # 6. summary
    log("\n" + "=" * 60)
    log(f"Source          : {label}")
    log(f"Runs resolved   : {len(runs)}")
    log(f"Files expected  : {n_files}")
    log(f"Files missing   : {len(missing)}")
    log(f"MD5 failures    : {'skipped' if opts.no_verify else len(failures)}")
    log(f"Mapping file    : {map_path}")
    log(f"Run report      : {report}")
    log(f"Finished        : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 60)

    for msg in failures:
        log("  " + msg)

    if missing or failures or opts.no_verify:
        log(f"{label}: INCOMPLETE - inspect the report and re-run "
            f"(existing files are skipped)")
        return 1

    marker = opts.outdir / f"{label}.download_complete"
    marker.write_text(
        f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{len(runs)} runs\t"
        f"{n_files} files\tMD5 verified\n"
    )
    log(f"{label}: SUCCESS - {n_files} file(s) downloaded and MD5-verified")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        die("interrupted", 130)