#!/usr/bin/env perl
# Creates a pass/fail file for easy to read results

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use List::Util qw/sum/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "3.0";
our $CITATION="SneakerNet pass/fail by Lee Katz";

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version tempdir=s citation check-dependencies help debug force numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
    }, $settings,
  );

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
  my @failHeader=sort keys(%{ $$failHash{$sample[0]} });

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
  print $failFh "#Failure by Kraken is when the number of contaminant reads is > 10%\n";
  close $failFh;

  return $failFile;
}

sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %whatFailed=();  # reasons why it failed

  # Read metrics
  my %readMetrics = ();
  # If the readMetrics.tsv file exists, read it.
  if(-e "$dir/readMetrics.tsv"){
    # Read the readMetrics file into %readMetrics
    open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
    my @header=split(/\t/,<READMETRICS>); chomp(@header);
    while(<READMETRICS>){
      chomp;
      my %F;
      @F{@header}=split(/\t/,$_);
      $F{File} = basename($F{File});
      $readMetrics{$F{File}} = \%F;
    }
    close READMETRICS;
  }
  # If the readmetrics file does not exist, then the values are blank
  else {
    logmsg "WARNING: readMetrics.tsv was not found. Unknown whether samples pass or fail by coverage and quality.";
  }

  # Kraken results
  my %kraken;
  if(-e "$dir/SneakerNet/forEmail/kraken.tsv"){
    open(my $krakenFh, '<', "$dir/SneakerNet/forEmail/kraken.tsv") or logmsg "WARNING: kraken results were not found in $dir/SneakerNet/forEmail/kraken.tsv: $!";
    my $header = <$krakenFh>;
    chomp($header);
    my @header = split(/\t/, $header);
    while(<$krakenFh>){
      chomp;
      next if(/^\s*#/);
      my @F = split(/\t/, $_);
      my %F;
      @F{@header} = @F;
      $F{PERCENTAGE_CONTAMINANT} //= 0;

      $kraken{$F{NAME}} = \%F;
    }
    close $krakenFh;
  }

  # Understand for each sample whether it passed or failed
  # on each category
  for my $samplename(keys(%$sampleInfo)){

    # Possible values for each:
    #  1: failed the category
    #  0: passed (did not fail)
    # -1: unknown
    my %fail=(
      coverage=>-1,
      quality => 0, # by default passes
      kraken  => -1,
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
      $$fastqMetrics{coverage} //= '.';    # by default, dot.
      if($$fastqMetrics{coverage} eq '.'){ # dot means coverage is unknown
        $totalCoverage = -1; # -1 means 'unknown' coverage
      } else {
        $$fastqMetrics{coverage} ||= 0; # force it to be a number if it isn't already
        $totalCoverage += $$fastqMetrics{coverage};
        logmsg "Sample $samplename += $$fastqMetrics{coverage}x => ${totalCoverage}x" if($$settings{debug});
      }

      $$sampleInfo{$samplename}{taxonRules}{quality} ||= 0;
      # If avgQual is missing, then -1 for unknown pass status
      $$fastqMetrics{avgQuality} //= '.';
      if($$fastqMetrics{avgQuality} eq '.'){
        $is_passing_quality{$fastq} = -1;
      }
      # Set whether this fastq passes quality by the > comparison:
      # if yes, then bool=true, if less than, bool=false
      else {
        $is_passing_quality{$fastq} =  $$fastqMetrics{avgQuality} >= $$sampleInfo{$samplename}{taxonRules}{quality};
      }
    }

    # Set whether the sample fails coverage
    #logmsg "DEBUG"; $$sampleInfo{$samplename}{taxonRules}{coverage} = 5;
    if($totalCoverage < $$sampleInfo{$samplename}{taxonRules}{coverage}){
      logmsg "  I will fail this sample $samplename" if($$settings{debug});
      $fail{coverage} = 1;
    } else {
      logmsg "  I will not fail this sample $samplename" if($$settings{debug});
      $fail{coverage} = 0;
    }
    if($totalCoverage == -1){
      logmsg "  ==> -1" if($$settings{debug});
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
        $fail{quality} = -1;
        last;
      }
    }

    # Did it pass by kraken / read classification
    if(defined($kraken{$samplename})){
      if($kraken{$samplename}{PERCENTAGE_CONTAMINANT} > 10){
        $fail{kraken} = 1;
      } else {
        $fail{kraken} = 0;
      }
    } else {
      $fail{kraken} = -1;
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

