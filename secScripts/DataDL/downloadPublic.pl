#!/usr/bin/env perl
#
# downloadPublic.pl - MATAFILER front end for retrieving public metagenomes.
#
# Resolves a study accession against the ENA portal, downloads the FASTQ files
# and writes a MATAFILER mapping file that can be fed straight into MATAF4.pl.
# The actual work is done by secScripts/DataDL/fetchPublicReads.py; this
# wrapper exists to take program paths from the MATAFILER config and to hand
# the job to the cluster in the way the rest of MATAFILER does.
#
# NETWORK, PLEASE READ:
#   This tool needs outbound HTTPS to www.ebi.ac.uk. On many clusters the
#   ordinary compute nodes have none, so a download submitted as a normal job
#   dies with a misleading "no runs resolved" error.
#
#   MATAFILER's job system (Mods::Subm) knows about short/medium/long/gpu/
#   highMem queues, but it has no notion of a queue that can reach the
#   internet. Until that gap is resolved in Subm.pm, this script offers two
#   supported routes:
#
#     default        run here, on the node you launched from. Login and
#                    submit nodes normally do have network. Fine for a
#                    dry run or a small study; wrap it in screen/tmux for
#                    anything large.
#
#     -sbatchOut     write a standalone sbatch script targeting a partition
#                    you name (-partition), which you then submit yourself.
#                    Use this for real downloads on clusters that provide a
#                    dedicated network-capable partition.
#
#   -submit 1 hands the job to MATAFILER's own queue system. It will only work
#   if your DEFAULT queue happens to have internet access. It usually does not.
#
# Usage:
#   ./downloadPublic.pl -acc PRJNA703757 -outDir /scratch/raw/PRJNA703757 \
#       -map /path/PRJNA703757.map -sampleTable PRJNA703757.sample_table.tsv \
#       -outPath /path/matafiler_output -relaxSmplID 1
#
#   Dry run first. It resolves the study and writes the mapping file WITHOUT
#   downloading anything, so you can inspect what you are about to fetch:
#   ./downloadPublic.pl -acc PRJNA703757 -outDir /scratch/x -map /tmp/x.map -dryRun 1
#
use warnings;
use strict;
use Getopt::Long qw( GetOptions );
use File::Basename;
use Cwd 'abs_path';

use Mods::IO_Tamoc_progs qw(getProgPaths setConfigFile);
use Mods::Subm qw(qsubSystem emptyQsubOpt);

sub usage {
    print STDERR <<"END";

downloadPublic.pl -- fetch public metagenomes from ENA/SRA and emit a
                     MATAFILER mapping file.

Required:
  -acc <ID>            study accession (PRJNA..., PRJEB..., SRP..., ERP...)
     or
  -metadata <file>     a pre-downloaded ENA filereport TSV instead of -acc

  -outDir <dir>        where the FASTQ files go
  -map <file>          mapping file to write

Study design (strongly recommended):
  -sampleTable <file>  TSV: run_accession, SmplID, AssmblGrps. Without it,
                       SmplID defaults to the run accession and nothing is
                       co-assembled. It also acts as a whitelist: only the
                       runs it lists are fetched.
  -outPath <dir>       written to the #OutPath header of the mapping file
  -relaxSmplID <0|1>   write #RelaxSMPLID TRUE            [0]

Run selection:
  -excludeStrategy <S> drop these ENA library_strategy values, comma
                       separated. AMPLICON is the usual one: many studies
                       deposit 16S and shotgun runs under one accession, and
                       assembling amplicon reads fails silently.
  -libStrategy <S>     keep only these library_strategy values. Careful:
                       the field is submitter-declared and often wrong.
  -minReads <N>        drop runs below N reads (blanks, failed libraries).
                       Runs declaring 0 reads are always dropped.  [0]

Execution:
  -cores <N>           parallel download streams. Keep low, ENA throttles. [4]
  -dryRun <0|1>        resolve and write the mapping file, download nothing.
                       ALWAYS do this first: a real run of a large study
                       moves hundreds of gigabytes.                     [0]

  -sbatchOut <file>    do not download; write a self-contained sbatch script
                       here, for you to submit. The recommended route for
                       real downloads on clusters with a network-capable
                       partition.
  -partition <name>    partition named in the sbatch script (e.g. a dedicated
                       download queue)
  -time <hh:mm:ss>     wall time in the sbatch script            [96:00:00]
  -mem <NG>            memory                                          [4G]

  -submit <0|1>        hand the job to MATAFILER's queue system. Works only
                       if your default queue has outbound HTTPS, which is
                       uncommon. Prefer -sbatchOut.                     [0]

END
    exit 1;
}

my ($acc, $metadata, $outDir, $mapOut, $sampleTable, $outPath) = ("", "", "", "", "", "");
my ($excludeStrategy, $libStrategy) = ("", "");
my ($sbatchOut, $partition) = ("", "");
my ($relaxSmplID, $minReads, $cores, $dryRun, $submit) = (0, 0, 4, 0, 0);
my $mem = "4G";
my $wallTime = "96:00:00";

GetOptions(
    "acc=s"              => \$acc,
    "metadata=s"         => \$metadata,
    "outDir=s"           => \$outDir,
    "map=s"              => \$mapOut,
    "sampleTable=s"      => \$sampleTable,
    "outPath=s"          => \$outPath,
    "relaxSmplID=i"      => \$relaxSmplID,
    "excludeStrategy=s"  => \$excludeStrategy,
    "libStrategy=s"      => \$libStrategy,
    "minReads=i"         => \$minReads,
    "cores=i"            => \$cores,
    "dryRun=i"           => \$dryRun,
    "submit=i"           => \$submit,
    "sbatchOut=s"        => \$sbatchOut,
    "partition=s"        => \$partition,
    "time=s"             => \$wallTime,
    "mem=s"              => \$mem,
) or usage();

usage() if (!$acc && !$metadata);
usage() if ($acc && $metadata);
usage() if (!$outDir || !$mapOut);

# getProgPaths dies with a clear message if the key is absent from the config,
# which is what we want: a missing entry is a broken install, not something to
# paper over with a guessed path.
my $fetchScr = getProgPaths("fetchReads_scr");
my $wgetBin  = getProgPaths("wget");

my $label = $acc ? $acc : basename($metadata);
$label =~ s/\.[^.]+$//;

my @args;
push @args, ($acc ? ("-a", $acc) : ("-m", $metadata));
push @args, ("-o", $outDir);
push @args, ("--map-out", $mapOut);
push @args, ("--run-id", $label);
push @args, ("--wget", $wgetBin);
push @args, ("-t", $cores);
push @args, ("--sample-table", $sampleTable)        if ($sampleTable);
push @args, ("--out-path", $outPath)                if ($outPath);
push @args, ("--relax-smplid")                      if ($relaxSmplID);
push @args, ("--exclude-strategy", $excludeStrategy) if ($excludeStrategy);
push @args, ("--library-strategy", $libStrategy)    if ($libStrategy);
push @args, ("--min-reads", $minReads)              if ($minReads > 0);
push @args, ("--dry-run")                           if ($dryRun);

my $cmd = "$fetchScr " . join(" ", @args) . "\n";

# --- write a standalone sbatch script -----------------------------------------
if ($sbatchOut) {
    my $partLine = $partition ? "#SBATCH --partition=$partition\n" : "";
    if (!$partition) {
        print STDERR "NOTE: no -partition given. The download will land on "
                   . "your cluster's default partition, which may have no "
                   . "internet access.\n";
    }
    open my $fh, ">", $sbatchOut or die "cannot write $sbatchOut: $!\n";
    print $fh <<"SB";
#!/bin/bash
#SBATCH --job-name=DL_$label
$partLine#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$cores
#SBATCH --mem=$mem
#SBATCH --time=$wallTime
#SBATCH --output=$outDir/DL_$label.%j.log
#SBATCH --error=$outDir/DL_$label.%j.err

# Generated by MATAFILER downloadPublic.pl
# Resumable: wget -c skips whatever is already on disk, so re-submitting
# after a timeout or a failure picks up where it stopped.

set -uo pipefail
mkdir -p $outDir

$cmd
SB
    close $fh;
    chmod 0755, $sbatchOut;
    print "sbatch script written: $sbatchOut\n";
    print "  partition : " . ($partition ? $partition : "<cluster default>") . "\n";
    print "  cores     : $cores\n  mem       : $mem\n  walltime  : $wallTime\n\n";
    print "Submit it with:\n  sbatch $sbatchOut\n";
    exit 0;
}

if (!$submit) {
    # Run here. The launching node is normally a login/submit node and can
    # reach ENA; compute nodes frequently cannot.
    print "Running locally (compute nodes usually have no outbound HTTPS).\n";
    print "Use -submit 1 -qpartition <net-capable partition> to submit instead.\n\n";
    print "$cmd\n";
    my $rc = system($cmd);
    if ($rc != 0) {
        die "\ndownloadPublic.pl: fetchPublicReads.py exited with status "
          . ($rc >> 8) . ".\nIf every accession failed to resolve, check that "
          . "this node can reach https://www.ebi.ac.uk .\n";
    }
    print "\nDone. Mapping file: $mapOut\n";
    print "Feed it to MATAF4.pl once you are happy with it.\n";
    exit 0;
}

# --- cluster submission ------------------------------------------------------
# Only sane where compute nodes have internet access.
my $logDir = "$outDir/LOGandSUB/";
system "mkdir -p $logDir" unless (-d $logDir);

my ($jobName) = ("DL_$label");
my $QSBoptHR = emptyQsubOpt(1, "");

# Mods::Subm has no concept of a network-capable queue, so this goes to the
# default (medium) queue. On most clusters that queue cannot reach the
# internet. -sbatchOut exists precisely because of this gap.
print "Submitting $jobName to MATAFILER's default queue.\n";
print "WARNING: MATAFILER cannot currently target a network-capable queue.\n";
print "         If the job fails with 'no runs resolved', that is a network\n";
print "         restriction, not a bad accession. Use -sbatchOut -partition\n";
print "         <net-capable partition> instead.\n\n";

# args: 1[bash/err/out stem] 2[cmd] 3[cores] 4[RAM per core] 5[job name]
#       6[job deps] 7[cwd?] 8[submit now?] 9[restricted hosts] 10[opts]
my ($jobd, $tmpCmd) = qsubSystem(
    $logDir . "downloadPublic.sh", $cmd, $cores, $mem, $jobName,
    "", "", 1, [], $QSBoptHR
);

print "Submitted. Mapping file will appear at: $mapOut\n";
exit 0;