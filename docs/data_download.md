# Public data download

MATAFILER can retrieve public metagenomes from ENA/SRA and produce a mapping
file for them, so that a published BioProject can go from an accession to an
assembled dataset without hand-editing anything.

```
secScripts/DataDL/
  downloadPublic.pl      Perl front end. Start here.
  fetchPublicReads.py    resolve, download, verify, write the mapping file
  map2sampleTable.py     turn an existing mapping file back into a sample table
  auditMapIdentity.py    check an existing mapping file against ENA
  describeRuns.py        ask ENA what a set of runs actually is
```

Only the Python standard library and `wget` are needed. No conda environment,
no new dependency.

---

## Quick start

**Always dry-run first.** A real run of a large study moves hundreds of
gigabytes. `-dryRun 1` resolves the study and writes the mapping file *without
downloading anything*, so you can look before you leap.

```bash
./secScripts/DataDL/downloadPublic.pl \
  -acc PRJEB26842 \
  -outDir /scratch/raw/PRJEB26842 \
  -map /path/to/mapfiles/PRJEB26842.map \
  -dryRun 1
```

Inspect the mapping file. When you are happy with it, drop `-dryRun`.

---

## The three traps

These are not hypothetical. All three were found in a 13-project meta-analysis
that had already been curated by hand. **None of them causes an error message.**
They are where this tool's real value lies.

### 1. One accession, two kinds of data

A study often deposits its 16S amplicon runs and its shotgun runs under the same
project accession.

| Project | Total runs | Amplicon | Shotgun |
|---|---|---|---|
| PRJEB26842 | 116 | 58 AMPLICON | 58 shotgun |
| PRJNA529487 | 708 | 674 AMPLICON | 34 shotgun |
| PRJNA743361 | 437 | 311 AMPLICON | 126 shotgun |
| SRP045211 | 217 | 108 AMPLICON | 109 shotgun |

Feed the amplicon runs into an assembler and it will not complain. It will
assemble them, and the result will be meaningless.

The tool prints the `library_strategy` breakdown before it does anything, and
warns loudly if a study mixes strategies and you asked for no selection.

### 2. `library_strategy` is whatever the submitter typed

It is not validated. Do not trust it as a whitelist:

- `PRJEB26842`'s shotgun runs are labelled `OTHER`
- `SRP045211`'s shotgun runs are labelled `WGA`

`--library-strategy WGS` on either of these silently selects **nothing**.

Excluding is robust where including is not: nobody mislabels 16S as anything
other than `AMPLICON`. So:

```bash
-excludeStrategy AMPLICON      # robust, works across studies
-libStrategy WGS               # fragile: check the breakdown first
```

Better still, pass a `-sampleTable` (below), which is authoritative and does not
depend on the submitter having filled in a field correctly.

### 3. ENA accumulates empty and broken runs

Projects are not frozen when they are published. `PRJNA703757` gained a run in
November 2025, four years after publication, carrying `read_count` 0 and no
files of any kind. It had found its way into a mapping file as a sample that
could never have any data.

Runs reporting zero reads are always dropped, with a warning. Runs that exist
but expose no downloadable files are reported and skipped. Use `-minReads` to
also exclude blanks and failed libraries.

---

## Study design: what the tool will not guess

A mapping file mixes two very different kinds of information.

**Derivable from ENA.** `SmplPrefix`, `#DirPath`, `#RunID`. The tool downloaded
the files; it knows their names.

**Not derivable from ENA.** `SmplID` and `AssmblGrps`. These encode experimental
design: what a sample is called biologically, and which runs are the same
subject and should therefore be co-assembled. ENA does not carry that reliably,
and guessing it from sample titles works silently and wrongly often enough to be
dangerous.

So the tool does not guess. Supply the design with `-sampleTable`:

```
run_accession	SmplID	AssmblGrps
ERR2587764	shotgunfmppost11	fmp11
ERR2587774	shotgunfmppre11	fmp11
ERR2587765	shotgunfmppost14	fmp14
```

Runs sharing an `AssmblGrps` value are co-assembled — here, the pre- and
post-intervention samples of the same subject.

The table also acts as a **whitelist**: only the runs it lists are fetched. That
solves the amplicon problem exactly, without depending on `library_strategy`.

Without a sample table, `SmplID` defaults to the run accession and every sample
becomes its own assembly group (no co-assembly). Perfectly usable, just not
aware of your study.

If you already have a curated mapping file, you do not have to build the table
by hand:

```bash
./secScripts/DataDL/map2sampleTable.py PRJEB26842.map -o PRJEB26842.sample_table.tsv
```

---

## How integrity is checked

A successful download is not the same as a correct one: exit code 0 does not
mean you got the file you asked for. Downloading is checked in two stages.

**During download (per file), three cheap checks** — none reads the whole file,
so bad files are rejected fast:

1. empty (`size == 0`)
2. not gzip (first two bytes are not the gzip magic `1f 8b` — catches HTML error
   pages saved as `.fastq.gz`)
3. wrong size (on-disk bytes ≠ ENA's declared `fastq_bytes` — catches truncated
   downloads)

**After everything is downloaded, one MD5 pass** — the only check that reads
full file contents, catching corruption that slips past the cheap checks. This
is not cheap (~30 min for 750 GB) but it is the final backstop, and in one real
download it caught six corrupt files that would otherwise have gone straight
into assembly.

### Recovery: backoff between retries

ENA's mirror nodes hiccup periodically, and failures come in windows lasting
minutes to tens of minutes, during which the same files repeatedly fail or come
back truncated. `wget`'s own retries burn out within seconds, all inside the same
window, so retrying immediately is nearly useless.

So the tool waits a growing amount of time between retry passes (2, 4, 6, 8
minutes by default), deliberately waiting for the failure window to pass before
re-downloading the incomplete files. This is what lets a large download run
unattended: without it, every project of a few hundred samples needs a human to
patch files by hand.

Downloads are resumable (`wget -c`), so a job that times out can simply be
resubmitted; it picks up where it stopped.

---

## Network access

**This tool needs outbound HTTPS to `www.ebi.ac.uk`.**

On most clusters, ordinary compute nodes have no internet. A download submitted
as a normal job will fail there with a misleading "no runs resolved" — which
looks like a bad accession but is a firewall.

Three routes, in order of preference:

**1. `-sbatchOut` (recommended for real downloads)**

Writes a standalone sbatch script targeting a partition you name, which you then
submit yourself. Use this when your cluster has a dedicated network-capable
partition.

```bash
./secScripts/DataDL/downloadPublic.pl \
  -acc PRJNA703757 -outDir /scratch/raw/PRJNA703757 \
  -map PRJNA703757.map -sampleTable PRJNA703757.sample_table.tsv \
  -sbatchOut dl_PRJNA703757.sh -partition <net-capable partition> \
  -time 96:00:00

sbatch dl_PRJNA703757.sh
```

Note that generating the script and submitting it may need to happen on
different nodes if your cluster separates internet access from the scheduler
(resolving the study needs the network; `sbatch` needs the scheduler).

**2. Run it locally (the default)**

Runs on whichever node you launched from. Login and submit nodes usually do have
network. Fine for a dry run or a small study. Wrap it in `screen` or `tmux` for
anything large, or an SSH drop will kill it.

**3. `-submit 1`**

Hands the job to MATAFILER's own queue system. This only works if your *default*
queue has internet access, which is uncommon — `Mods::Subm` currently has no way
to express "this job needs the internet". Prefer `-sbatchOut`.

---

## Auditing what you already have

`auditMapIdentity.py` compares the run accessions in existing mapping files
against what ENA actually holds, per `library_strategy`.

```bash
./secScripts/DataDL/auditMapIdentity.py /path/to/mapfiles/ --out audit.tsv
```

For each project it reports:

- **in the map, not in ENA** — never benign. Typos, duplicates, withdrawn runs,
  or samples borrowed from another study.
- **in ENA, not in the map** — broken down by strategy. `AMPLICON` here is
  expected. Shotgun here needs an explanation.
- duplicate accessions across map rows.

Counting is not enough: "58 in the map, 58 shotgun in ENA" does not prove they
are the *same* 58. This compares the accession sets.

When it flags something, ask ENA what those runs actually are:

```bash
./secScripts/DataDL/describeRuns.py PRJNA675301 --not-in-map PRJNA675301.map
```

which is how the 40 runs missing from PRJNA675301 turned out to be extraction
blanks, empty-well controls and Vibrio/lambda spike-ins — correctly excluded,
but nobody had written that down.

---

## Options (`downloadPublic.pl`)

**Required**

| Flag | Meaning |
|---|---|
| `-acc <ID>` | study accession (PRJNA/PRJEB/SRP/ERP…) |
| `-metadata <file>` | …or a pre-downloaded ENA filereport TSV instead |
| `-outDir <dir>` | where the FASTQ files go |
| `-map <file>` | mapping file to write |

**Study design**

| Flag | Meaning |
|---|---|
| `-sampleTable <file>` | `run_accession` / `SmplID` / `AssmblGrps`. Also a whitelist. |
| `-outPath <dir>` | written to the `#OutPath` header |
| `-relaxSmplID <0\|1>` | write `#RelaxSMPLID TRUE` |

**Run selection**

| Flag | Meaning |
|---|---|
| `-excludeStrategy <S>` | drop these `library_strategy` values (e.g. AMPLICON) |
| `-libStrategy <S>` | keep only these — fragile, see trap 2 |
| `-minReads <N>` | drop runs below N reads; zero-read runs always go |

**Execution**

| Flag | Meaning |
|---|---|
| `-cores <N>` | parallel download streams. ENA throttles; keep low. Default 4 |
| `-dryRun <0\|1>` | resolve and write the map, download nothing |
| `-sbatchOut <file>` | write an sbatch script instead of downloading |
| `-partition <name>` | partition for that script |
| `-time <hh:mm:ss>` | wall time for that script. Default 96:00:00 |

The underlying `fetchPublicReads.py` additionally accepts `--resume-passes <N>`
(retry passes after the first, default 4) and `--resume-wait <sec>` (base wait
before pass N, actual wait = `resume-wait × N`, capped at 900s, default 120).

---

## Recommended workflow

```bash
MD=/path/to/project

# 1. what is actually in this study?
./secScripts/DataDL/downloadPublic.pl -acc PRJEB26842 \
    -outDir /scratch/raw/PRJEB26842 -map /tmp/preview.map -dryRun 1

# 2. build the sample table encoding your design (or convert an existing map)
./secScripts/DataDL/map2sampleTable.py old.map -o $MD/PRJEB26842.sample_table.tsv

# 3. dry-run with the design, and read the mapping file it produces
./secScripts/DataDL/downloadPublic.pl -acc PRJEB26842 \
    -outDir /scratch/raw/PRJEB26842 -map $MD/mapfiles/PRJEB26842.map \
    -sampleTable $MD/PRJEB26842.sample_table.tsv \
    -outPath $MD/matafiler_output -relaxSmplID 1 \
    -excludeStrategy AMPLICON -dryRun 1

# 4. download for real
./secScripts/DataDL/downloadPublic.pl <same args> \
    -sbatchOut $MD/dl.sh -partition <net-capable partition>
sbatch $MD/dl.sh

# 5. hand the mapping file to MATAFILER
```
