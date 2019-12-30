#!/usr/bin/env perl
# Transfers files to a remote destination and QCs them beforehand.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use List::Util qw/sum/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "2.0";
our $CITATION="SneakerNet pass/fail by Lee Katz";

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version tempdir=s citation help debug force numcpus=i)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }
  if($$settings{citation}){
    print $CITATION."\n";
    return 0;
  }

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";

  my $outfile=passfail($dir,$settings);
  logmsg "The pass/fail file is under $outfile";
  
  recordProperties($dir,{version=>$VERSION,table=>$outfile});

  return 0;
}

sub passfail{
  my($dir,$settings)=@_;

  my $failFile="$dir/SneakerNet/forEmail/passfail.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my $failHash=identifyBadRuns($dir,$sampleInfo,$settings);
  
  my @sample=keys(%$failHash);
  my @failHeader=keys(%{ $$failHash{$sample[0]} });

  open(my $failFh, ">", $failFile) or die "ERROR: could not write to $failFile: $!";
  print $failFh join("\t", "Sample", @failHeader)."\n";
  while(my($sample,$fail)=each(%$failHash)){
    print $failFh $sample;
    for(@failHeader){
      print $failFh "\t".$$fail{$_};
    }
    print $failFh "\n";
  }
  print $failFh "#1: fail\n#0: pass\n#-1: unknown\n";
  close $failFh;

  return $failFile;
}

sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %whatFailed=();  # reasons why it failed

  # Read the readMetrics file into %readMetrics
  open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
  my @header=split(/\t/,<READMETRICS>); chomp(@header);
  my %readMetrics = ();
  while(<READMETRICS>){
    chomp;
    my %F;
    @F{@header}=split(/\t/,$_);
    $F{File} = basename($F{File});
    $readMetrics{$F{File}} = \%F;
  }
  close READMETRICS;

  # Understand for each sample whether it passed or failed
  # on each category
  for my $samplename(keys(%$sampleInfo)){

    # Possible values for each:
    #  1: failed the category
    #  0: passed (did not fail)
    # -1: unknown
    my %fail=(
      coverage=>-1,
      quality => 1, # by default passes
    );

    # Skip anything that says undetermined.
    if($samplename=~/^Undetermined/i){
      next;
    }

    # Get the name of all files linked to this file through the sample.
    $$sampleInfo{$samplename}{fastq}//=[]; # Set {fastq} to an empty list if it does not exist
    my @file=@{$$sampleInfo{$samplename}{fastq}};

    # Get metrics from the fastq files
    my $totalCoverage = 0;
    my %is_passing_quality = (); # Whether a read passes quality
    for my $fastq(@file){
      my $fastqMetrics = $readMetrics{basename($fastq)};

      # Coverage
      if($$fastqMetrics{coverage} eq '.'){ # dot means coverage is unknown
        $totalCoverage = -1; # -1 means 'unknown' coverage
      } else {
        $$fastqMetrics{coverage} ||= 0; # force it to be a number if it isn't already
        $totalCoverage += $$fastqMetrics{coverage};
      }

      # Set whether this fastq passes quality by the > comparison:
      # if yes, then bool=true, if less than, bool=false
      $$sampleInfo{$samplename}{taxonRules}{quality} ||= 0;
      $is_passing_quality{$fastq} =  $$fastqMetrics{avgQuality} >= $$sampleInfo{$samplename}{taxonRules}{quality};
    }

    # Set whether the sample fails coverage
    if($totalCoverage < $$sampleInfo{$samplename}{taxonRules}{coverage}){
      $fail{coverage} = 1;
    }
    if($totalCoverage == -1){
      $fail{coverage} = -1;
    }

    # if any filename fails quality, then overall quality fails
    while(my($fastq,$passed) = each(%is_passing_quality)){
      # If
      if($passed == 0){
        $fail{$samplename} = 1;
        last;
      }
      # Likewise if a file's quality is unknown, then it's all unknown
      if($passed == -1){
        $fail{$samplename} = -1;
        last;
      }
    }

    $whatFailed{$samplename} = \%fail;
  }

  return \%whatFailed;
}

sub usage{
  print "Passes or fails a run based on available information
  Usage: $0 MiSeq_run_dir
  --version
";
  exit(0);
}

