#!/usr/bin/env perl
# Reformats a sample sheet

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use Cwd qw/realpath/;
use File::Temp;
use File::Copy qw/cp/;
use FindBin;
use Config::Simple;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig samplesheetInfo passfail command logmsg version/;

my $snVersion=version();

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help numcpus=i debug tempdir=s force)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  # We have to know if the argument is a directory or file, and so
  # we need to resolve the link.
  if(-l $ARGV[0]){
    $ARGV[0] = readlink($ARGV[0]);
  }

  my $outfile     = "";
  my $samplesheet = "";
  if(-d $ARGV[0]){
    $outfile = "$ARGV[0]/samples.tsv";
    $samplesheet = findSampleSheet($ARGV[0], $settings);
  } else {
    $samplesheet = $ARGV[0];
  }

  my $sampleHash = samplesheetInfo($samplesheet,$settings);
  my $tsv        = sampleHashToTsv($sampleHash, $settings);

  # Print output to stdout and to output file
  print $tsv;
  if($outfile){
    if(-e $outfile && !$$settings{force}){
      logmsg "NOTE: will not write to $outfile: it already exists.  Overwrite with --force.";
    } else {
      open(my $outfh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
      print $outfh $tsv;
      close $outfh;
      logmsg "Wrote samples file to $outfile";
    }
  }
  
  return 0;
}

sub findSampleSheet{
  my($dir, $settings)=@_;

  if(-e "$dir/SampleSheet.csv"){
    return "$dir/SampleSheet.csv";
  }

  if(-e "$dir/SampleSheetUsed.csv"){
    cp("$dir/SampleSheetUsed.csv", "$dir/SampleSheet.csv");
    _removeRunNumberFromSamples("$dir/SampleSheet.csv", $settings);
    return "$dir/SampleSheet.csv";
  }

  die "ERROR: could not find the sample sheet file in $dir/";
}

# Edit a sample sheet in-place to remove a run identifier
# from the sample names. For some reason the Miniseq
# appends a four digit number, e.g. "-6006" to the end
# of each sample name.
sub _removeRunNumberFromSamples{
  my($samplesheet,$settings)=@_;

  my $newSamplesheetString="";
  open(SAMPLESHEET,"<", $samplesheet) or die "ERROR: could not read $samplesheet: $!";
  my $reachedSamples=0;
  my $runid="";
  while(<SAMPLESHEET>){
    # Make a note of the run ID when I see it
    if(/Local Run Manager Analysis Id,\s*(\d+)/){
      $runid=$1;
    }

    if(!$reachedSamples){
      $newSamplesheetString.=$_;
      if(/Sample_ID,/){
        $reachedSamples=1;
      }
    }
    # Read the samples and remove the run ID
    else {
      my($samplename,@therest)=split(/,/,$_);
      $samplename=~s/\-$runid$//;
      $newSamplesheetString.=join(",",$samplename,@therest);
    }
  }
  close SAMPLESHEET;

  # Now rewrite the sample sheet
  open(SAMPLESHEET,">", $samplesheet) or die "ERROR: could not write to $samplesheet: $!";
  print SAMPLESHEET $newSamplesheetString;
  close SAMPLESHEET;

  return 1;
}

sub sampleHashToTsv{
  my($sampleHash,$settings)=@_;

  my $tsv = "";
  for my $sample(sort{ $a cmp $b } keys(%$sampleHash)){
    next if(ref($$sampleHash{$sample}) ne 'HASH');
    my $fastq = $$sampleHash{$sample}{fastq};
    my $taxon = $$sampleHash{$sample}{taxonRules}{taxon};

    $tsv .= join("\t", 
      $sample,         # column1: sample name 
      $taxon,          # 2: describe what rules to follow
      @$fastq,         # 3-: fastq files in order of R1, R2
    )."\n";
  }
  
  return $tsv;
}

sub usage{
  "Reformat a SampleSheet.csv file in a SneakerNet run
  Usage: $0 run-dir|SampleSheet.csv
  
  Because this is a SneakerNet plugin, the first positional
  argument must be the SneakerNet run.  However, this script
  will also accept the spreadsheet itself as the first
  argument.

  OUTPUT: this script will print the new samplesheet to 
  stdout but if the first argument is a directory, it
  will place the file into the directory.
  "
}
