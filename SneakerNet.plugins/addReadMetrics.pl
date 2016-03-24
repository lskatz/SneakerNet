#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use Email::Stuffer;
use List::MoreUtils qw/uniq/;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig logmsg samplesheetInfo command/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug test numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  
  my $dir=$ARGV[0];

  addReadMetrics($dir,$settings);

  # Mark this file as something to attach for an email later
  symlink("$dir/readMetrics.tsv","$dir/SneakerNet/forEmail/readMetrics.tsv");

  return 0;
}

sub addReadMetrics{
  my($dir,$settings)=@_;

  return if(-e "$dir/readMetrics.tsv");

  logmsg "Reading sample $dir/SampleSheet.csv";
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  logmsg "Running fast read metrics";
  command("run_assembly_readMetrics.pl --numcpus $$settings{numcpus} --fast $dir/*.fastq.gz | sort -k3,3n > $dir/readMetrics.tsv.tmp");


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

# Use a hash of a line from a readmetrics output 
# to determine genome coverage.
sub calculateCoverage{
  my($h,$sampleInfo,$settings)=@_;

  # $h contains read metrics for this one row.

  my $file=basename($$h{File});
  my $samplename=$file || "";
  $samplename=~s/_S\d+_.*?$//;

  # Find out if this file has an expected genome size from the Sample Sheet.
  my $expectedGenomeSize=0;
  my $organism="";
  if($$sampleInfo{$samplename}{expectedgenomesize}){
    $expectedGenomeSize=$$sampleInfo{$samplename}{expectedgenomesize} * 10**6;
  }elsif($$sampleInfo{$samplename}{taxonRules}{genomeSize}){
    $expectedGenomeSize=$$sampleInfo{$samplename}{taxonRules}{genomeSize};
  }

  my $coverage=$$h{coverage} || 0; 

  # Recalculate coverage, if it's possible
  if($expectedGenomeSize > 0){
    $coverage=$$h{totalBases}/$expectedGenomeSize;
    $coverage=sprintf("%0.2f",$coverage); # round it
    logmsg "Decided that $$h{File} is $organism with expected genome size $expectedGenomeSize. Calculated coverage: $coverage";
  } else {
    logmsg "Warning: could not understand what organism $$h{File} belongs to; coverage was not recalculated. Reported coverage: $coverage";
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
  "
}

