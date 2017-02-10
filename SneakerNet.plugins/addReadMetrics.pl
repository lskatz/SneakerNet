#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

use threads;
use Thread::Queue;

use lib "$FindBin::RealBin/../lib/perl5";
use Email::Stuffer;
use List::MoreUtils qw/uniq/;
use SneakerNet qw/readConfig logmsg samplesheetInfo command/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help force inbox=s debug test numcpus=i tempdir=s)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir($0.".XXXXXX", TMPDIR=>1, CLEANUP=>1);
  logmsg "Tempdir is $$settings{tempdir}";
  
  my $dir=$ARGV[0];

  addReadMetrics($dir,$settings);

  # Mark this file as something to attach for an email later
  link("$dir/readMetrics.tsv","$dir/SneakerNet/forEmail/readMetrics.tsv");

  return 0;
}

sub addReadMetrics{
  my($dir,$settings)=@_;

  return if(-e "$dir/readMetrics.tsv" && !$$settings{force});

  logmsg "Reading sample $dir/SampleSheet.csv";
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  logmsg "Running fast read metrics";

  my $Q=Thread::Queue->new(glob("$dir/*.fastq.gz"));
  my @thr;
  for (0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&readMetricsWorker, $Q, $settings);
    $Q->enqueue(undef);
  }

  for(@thr){
    $_->join;
  }

  system("cat $$settings{tempdir}/*/readMetrics.tsv | head -n 1 > $dir/readMetrics.tsv.tmp"); # header
  system("sort -k3,3n $$settings{tempdir}/*/readMetrics.tsv | uniq -u >> $dir/readMetrics.tsv.tmp"); # content
  die if $?;

  # edit read metrics to include genome sizes
  logmsg "Backfilling values in $dir/readMetrics.tsv";
  my $newReadMetrics;
  open(READMETRICS,"$dir/readMetrics.tsv.tmp") or die "ERROR: could not open $dir/readMetrics.tsv.tmp because $!";
  open(READMETRICSFINAL,">","$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv for writing: $!";

  # get the header and also put it into the final output file
  my $header=<READMETRICS>;
  print READMETRICSFINAL $header;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<READMETRICS>){
    chomp;
    # read in each line into the appropriate header
    my %h;
    @h{@header}=split(/\t/);

    # find the genome size based on the filename
    my $coverage=calculateCoverage(\%h,$sampleInfo,$settings);
    $h{coverage}=$coverage;
    $h{File}=basename($h{File});

    for(@header){
      print READMETRICSFINAL "$h{$_}\t";
    }
    print READMETRICSFINAL "\n";
  }
  close READMETRICSFINAL;
  close READMETRICS;

  # Clean up by removing the temporary file
  unlink("$dir/readMetrics.tsv.tmp");
}

sub readMetricsWorker{
  my($Q, $settings)=@_;

  my $tempdir=tempdir("worker.XXXXXX", DIR=>$$settings{tempdir}, CLEANUP=>1);
  while(defined(my $fastq=$Q->dequeue)){
    logmsg "read metrics for $fastq";
    eval{
      command("run_assembly_readMetrics.pl --numcpus 1 --fast $fastq >> $tempdir/readMetrics.tsv");
      return 1;
    };
    if($@){
      logmsg "There was an error running run_assembly_readMetrics.pl.  This might be because of a divide-by-zero error. This can be solved by running the metrics without subsampling the reads which is slower.\n";
      logmsg "Rerunning without --fast.";
      command("run_assembly_readMetrics.pl --numcpus 1 $fastq >> $tempdir/readMetrics.tsv");
    } 

  }
}


# Use a hash of a line from a readmetrics output 
# to determine genome coverage.
sub calculateCoverage{
  my($h,$sampleInfo,$settings)=@_;

  # $h contains read metrics for this one row.

  my $file=basename($$h{File});
  my $samplename=$file || "";
  $samplename=~s/_S\d+_.*?$//;
  # Parse HiSeq names correctly
  if($file=~/(^.+_SAN\d+\w\d+)/){ # e.g., 2015V-1036_SAN5927A15_TCCGGAGA-CAGGACGT_L001_R1_001.fastq.gz
    $samplename=$1;               #       2015V-1036_SAN5927A15
  }

  # Find out if this file has an expected genome size from the Sample Sheet.
  my $expectedGenomeSize=0;
  my $organism="";
  if($$sampleInfo{$samplename}{expectedgenomesize}){
    $expectedGenomeSize=$$sampleInfo{$samplename}{expectedgenomesize} * 10**6;
  }elsif($$sampleInfo{$samplename}{taxonRules}{genomesize}){
    $expectedGenomeSize=$$sampleInfo{$samplename}{taxonRules}{genomesize};
  }

  my $coverage=$$h{coverage} || 0; 

  # Recalculate coverage, if it's possible
  if($expectedGenomeSize > 0){
    $coverage=$$h{totalBases}/$expectedGenomeSize;
    $coverage=sprintf("%0.2f",$coverage); # round it
    logmsg "Decided that $$h{File} is $organism with expected genome size $expectedGenomeSize. Calculated coverage: $coverage";
  } else {
    logmsg "Warning: could not understand what organism $$h{File} belongs to. I tried to look it up by $samplename. Coverage was not recalculated.";
  }
  return $coverage;
}
 
################
# Utility subs #
################

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 runDir
  --debug # Show debugging information
  --numcpus  1
  --tempdir ''
  --force
  "
}

