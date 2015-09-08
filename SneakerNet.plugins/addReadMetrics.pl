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

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug test)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  
  my $dir=$ARGV[0];

  addReadMetrics($dir,$settings);

  return 0;
}

sub addReadMetrics{
  my($dir,$settings)=@_;
  logmsg "Running fast read metrics";
  command("run_assembly_readMetrics.pl --fast $dir/*.fastq.gz | sort -k3,3n > $dir/readMetrics.tsv.tmp");

  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  # edit read metrics to include genome sizes
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

  # Two variables to figure out for genome coverage
  my $expectedGenomeSize=0;
  my $organism="";

  # Find out if this file has an expected genome size from the Sample Sheet.
  # This isn't the most efficient way to do it but then again this hash isn't too big.
  my $file=basename($$h{File});
  while(my($name,$info)=each(%$sampleInfo)){
    if(grep(/\Q$file\E/,@{ $$info{fastq} })){
      $expectedGenomeSize=(keys(%{ $$info{expectedgenomesize} }))[0] * 1e6;
      $organism=(keys(%{ $$info{species} }))[0];
    }
  }

  my $coverage=$$h{coverage} || '.'; # default value in case one isn't found

  # See if we can recalculate the coverage based on the filename
  if(!$expectedGenomeSize){
    # See if this filename matches any organism regex
    for my $info(@{ $$settings{genomeSizes} }){
      my($regex,$expectedGenomeSizeFromConfig);
      ($regex,$expectedGenomeSizeFromConfig,$organism)=@$info;

      # Calculate coverage from either the config file, the SampleSheet,
      # or else you just can't calculate.
      if($file=~/$regex/){
        $expectedGenomeSize=$expectedGenomeSizeFromConfig;
        last;
      }
    }
  }

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
 
sub samplesheetInfo{
  my($samplesheet,$settings)=@_;

  my $section="";
  my @header=();
  my %sample;
  open(SAMPLE,$samplesheet) or die "ERROR: could not open sample spreadsheet $samplesheet: $!";
  while(<SAMPLE>){
    chomp;
    if(/^\[(\w+)\]/){
      $section=lc($1);
      my $header=<SAMPLE>;
      chomp($header);
      @header=split(/,/,lc($header));
      next;
    }
    if($section eq "data"){
      my %F;
      @F{@header}=split(/,/,$_);
      for my $keyvalue(split(/;/,lc($F{description}))){
        my($key,$value)=split(/=/,$keyvalue);
        $key=~s/^\s+|\s+$//g;      #whitespace trim
        $value=~s/^\s+|\s+$//g;    #whitespace trim
        $F{$key}={} if(!$F{$key});
        $F{$key}{$value}++;
      }
      
      $sample{$F{sample_id}}=\%F;
    }
  }

  # Try to associate samples to files
  while(my($samplename,$sampleinfo)=each(%sample)){
    my @possibleFastq=glob(dirname($samplesheet)."/$samplename*.fastq.gz");
    $sample{$samplename}{fastq}=\@possibleFastq;
  }

  return \%sample;
}

################
# Utility subs #
################
sub readConfig{
  my @file=glob("$FindBin::RealBin/config/*");
  my $settings={};
  for(@file){
    open(CONFIGFILE,$_) or die "ERROR: could not open config file $_: $!";
    my $key=basename $_;
    while(<CONFIGFILE>){
      s/^\s+|\s+$//g; # trim
      next if(/^$/);
      next if(/^#/);
      my $configLine=[split(/\t/,$_)];
      push(@{ $$settings{$key} },$configLine);
    }
    close CONFIGFILE;
  }
  return $settings;
}


sub command{
  my($command,$settings)=@_;
  logmsg "COMMAND\n  $command" if($$settings{debug});
  system($command);
  die "ERROR running command\n  $command" if $?;
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 runDir
  --debug # Show debugging information
  "
}
