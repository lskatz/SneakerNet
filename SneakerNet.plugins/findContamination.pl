#!/usr/bin/env perl
# Use Kraken to find contamination.
# Modeling the workflow from https://github.com/lskatz/lskScripts/blob/205c03d3c4c44bfa2f961015c8d5f6971e5698d2/qsub/launch_kraken.sh

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;

my $KRAKENDIR="/opt/kraken";
my $KRONADIR="/opt/KronaTools-2.6/bin";
$ENV{PATH}="$ENV{PATH}:$KRAKENDIR:$KRONADIR";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help tempdir=s numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{KRAKEN_DEFAULT_DB} ||= die "ERROR: KRAKEN_DEFAULT_DB needs to be defined under config/settings";
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];
  
  my $outdir=runKrakenOnDir($dir,$settings);

  return 0;
}

sub runKrakenOnDir{
  my($dir,$settings)=@_;
  my $outdir="$dir/kraken";
  mkdir $outdir;

  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    my $sampledir="$outdir/$sampleName";
    mkdir $sampledir;
    runKraken($s,$sampledir,$settings);

    my $expectedSpecies=$$s{species} || "";
    my $percentContaminated=reportContamination($sampledir,$expectedSpecies,$settings);

    if($percentContaminated > 10){
      logmsg "$sampleName (taxon: $expectedSpecies) is $percentContaminated% contaminated";
    }
  }
}

sub runKraken{
  my($sample,$sampledir,$settings)=@_;

  my $html="$sampledir/report.html";
  return if(-e $html);

  my $reads=join(" ",@{ $$sample{fastq} });
  
  command("$KRAKENDIR/kraken --fastq-input --paired --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out $reads");

  command("$KRAKENDIR/kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  # TODO kraken-report

  command("$KRONADIR/ktImportText -o $html $sampledir/kraken.taxonomy");
}

sub reportContamination{
  my($sampledir,$expectedSpecies,$settings)=@_;
  return 100 if(!$expectedSpecies);

  my $taxfile="$sampledir/kraken.taxonomy";

  my $numCorrectReads=0;
  my $numContaminantReads=0;
  open(TAXONOMY,'<',$taxfile) or die "ERROR: could not open $taxfile for reading: $!";
  while(<TAXONOMY>){
    chomp;
    my($numReads,undef, $domain, $kingdom, $phylum, $class, $order, $family, $genus, $species)=split /\t/;
    $genus||="";
    $species||="";
    my $scientificName="$genus $species";

    if($expectedSpecies eq $genus || $expectedSpecies eq $species || $expectedSpecies eq $scientificName){
      $numCorrectReads+=$numReads;
    } else {
      $numContaminantReads+=$numReads;
    }
  }
  
  my $percentContamination=$numContaminantReads/($numContaminantReads+$numCorrectReads) * 100;
  return $percentContamination;
}

sub usage{
  "Finds contamination in a miseq run
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  "
}

