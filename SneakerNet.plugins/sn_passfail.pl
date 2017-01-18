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
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help debug force numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  my $outfile=passfail($dir,$settings);
  logmsg "The pass/fail file is under $outfile";
  
  return 0;
}

sub passfail{
  my($dir,$settings)=@_;

  my $failFile="$dir/SneakerNet/forEmail/passfail.tsv";
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my $failHash=identifyBadRuns($dir,$sampleInfo,$settings);
  
  my @sample=keys(%$failHash);
  my @failHeader=keys($$failHash{$sample[0]});

  open(my $failFh, ">", $failFile) or die "ERROR: could not write to $failFile: $!";
  print $failFh join("\t", "Sample", @failHeader)."\n";
  while(my($sample,$fail)=each(%$failHash)){
    print $failFh $sample;
    for(@failHeader){
      print $failFh "\t".$$fail{$_};
    }
    print $failFh "\n";
  }
  print $failFh "1: fail\n0: pass\n-1: unknown\n";
  close $failFh;

  return $failFile;
}

sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %whatFailed=();  # reasons why it failed

  open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
  my @header=split(/\t/,<READMETRICS>); chomp(@header);
  while(<READMETRICS>){
    chomp;
    my %F;
    @F{@header}=split(/\t/,$_);
    my $samplename=basename($F{File},'.fastq.gz');
    $samplename=~s/_S\d+_.*//; # figure out the sample name before the the _S1_ pattern

    # Possible values for each:
    #  1: failed the category
    #  0: passed (did not fail)
    # -1: unknown
    my %fail=(
      coverage=>-1,
      quality => 0,
    );

    # Skip anything that says undetermined.
    if($samplename=~/^Undetermined/){
      next;
    }

    # Get the name of all files linked to this file through the sample.
    $$sampleInfo{$samplename}{fastq}//=[];
    my @file=@{$$sampleInfo{$samplename}{fastq}};
    
    # Compare coverage of one read against half of the
    # threshold coverage because of PE reads.
    if($F{coverage} ne '.'){ # dot means coverage is unknown.
      $F{coverage}||=0;
      $$sampleInfo{$samplename}{taxonRules}{coverage}||=0;
      if($F{coverage} < $$sampleInfo{$samplename}{taxonRules}{coverage}/2){
        $fail{coverage}=1;
        logmsg "Low coverage in $F{File}";
      } else {
        $fail{coverage}=0;
      }
    }

    $F{avgQuality} ||= 0;
    $$sampleInfo{$samplename}{taxonRules}{quality}||=0;
    if($F{avgQuality} < $$sampleInfo{$samplename}{taxonRules}{quality}){
      $fail{quality}=1;
      logmsg "low quality in $F{File}";
    }
    $whatFailed{$F{File}}=\%fail;

  }

  return \%whatFailed;
}

sub usage{
  "Passes or fails a run based on available information
  Usage: $0 MiSeq_run_dir
  "
}

