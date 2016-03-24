#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  transferFilesToRemoteComputers($dir,$settings);

  return 0;
}

sub transferFilesToRemoteComputers{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  # Which files should be skipped according to Q/C?
  my $toSkip=identifyBadRuns($dir,$sampleInfo,$settings);

  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases
    my $taxon=$$s{species} || 'NOT LISTED';
    logmsg "The taxon of $sampleName is $taxon";
    if(grep {/calcengine/i} @{ $$s{route} }){
      for(@{ $$s{fastq} }){
        next if($$toSkip{basename($_)});
        my $subfolder=$$s{taxonRules}{dest_subfolder} || "SneakerNet";
        $filesToTransfer{$subfolder}.=$_." ";
      }
      logmsg "One route for sample $sampleName is the Calculation Engine";
    } else {
      logmsg "Note: The route for $sampleName was not listed in the sample sheet.";
    }
  }

  #die "ERROR: no files to transfer" if (!$filesToTransfer);
  logmsg "WARNING: no files will be transferred" if(!keys(%filesToTransfer));

  # Make the transfers based on taxon.
  while(my($subfolder,$fileString)=each(%filesToTransfer)){

    logmsg "Transferring to $subfolder:\n  $fileString";
    command("rsync --update -av $fileString $$settings{transfer_destination_string}/$subfolder/");
  }
}

# Which runs should be skipped based on bad quality
# or bad coverage, or whatever?
# Returns filenames.
sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %toSkip=();

  open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
  my @header=split(/\t/,<READMETRICS>); chomp(@header);
  while(<READMETRICS>){
    chomp;
    my %F;
    @F{@header}=split(/\t/,$_);
    my $samplename=basename($F{File},'.fastq.gz');
    $samplename=~s/_S\d+_.*//; # figure out the sample name before the the _S1_ pattern

    # Skip anything that says undetermined.
    if($samplename=~/^Undetermined/){
      $toSkip{basename($F{File})}=1;
      next;
    }

    #die Dumper [$samplename,$$sampleInfo{$samplename},\%F];

    # Compare coverage of one read against half of the
    # threshold coverage because of PE reads.
    if($F{coverage} ne '.'){ # dot means coverage is unknown.
      if($F{coverage} < $$sampleInfo{$samplename}{taxonRules}{coverage}/2){
        $toSkip{basename($F{File})}=1;
        logmsg "Low coverage in $F{File}";
      }
    }

    if($F{avgQuality} < $$sampleInfo{$samplename}{taxonRules}{quality}){
      my @file=@{$$sampleInfo{$samplename}{fastq}};
      $toSkip{basename($_)}=1 for(@file);
      logmsg "low quality in $F{File}\n  Skipping @file";
    }

  }
  
  return \%toSkip;
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 MiSeq_run_dir
  "
}

