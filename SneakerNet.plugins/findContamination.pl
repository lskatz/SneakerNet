#!/usr/bin/env perl
# Use Kraken to find contamination.
# Modeling the workflow from https://github.com/lskatz/lskScripts/blob/205c03d3c4c44bfa2f961015c8d5f6971e5698d2/qsub/launch_kraken.sh

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
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

  # make the report emailable 
  symlink(rel2abs("$outdir/report.tsv"), "$dir/SneakerNet/forEmail/kraken.tsv");

  return 0;
}

sub runKrakenOnDir{
  my($dir,$settings)=@_;
  my $outdir="$dir/SneakerNet/kraken";
  system("mkdir -p $outdir");

  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my %filesToTransfer=(); # hash keys are species names
  my @report; # reporting contamination in an array, in case I want to sort it later
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    # Skip any samples without reads, ie, samples that are misnamed or not sequenced.
    # There is no way to predict how a sample is misnamed and so it does not fall under
    # this script's purview.
    if(!defined($$s{fastq}) || !@{ $$s{fastq} }){
      logmsg "WARNING: I could not find the reads for $sampleName. Skipping.";
      next;
    }

    my $sampledir="$outdir/$sampleName";
    system("mkdir -p $sampledir");
    logmsg "Running Kraken on $sampleName";
    runKraken($s,$sampledir,$settings);

    my $expectedSpecies=$$s{species} || "";
    my ($percentContaminated,$html,$bestGuess)=reportContamination($sampledir,$expectedSpecies,$settings);
    next if(!$html);

    # Add onto the contamination report
    push(@report,join("\t",$sampleName,$expectedSpecies,$bestGuess,$percentContaminated));
    symlink(rel2abs($html),"$dir/SneakerNet/forEmail/$sampleName.kraken.html");

    # Report anything with >10% contamination to the printout.
    #if($percentContaminated > 10){
      #logmsg "$sampleName (taxon: $expectedSpecies) is $percentContaminated% contaminated";
    #}
  }

  # print the report to a file.
  # Note: Taylor reports 20% unclassified and a 4% contamination in one genome, which was considered uncontaminated.  Further testing is needed.
  unshift(@report,join("\t",qw(NAME LABELED_TAXON BEST_GUESS PERCENTAGE_OF_GENOME_IS_BEST_GUESS)));
  push(@report,"NOTE: A genome with >90% reads in agreement with its species has not been shown to indicate contamination");
  push(@report,"NOTE: More testing needs to be performed before any conclusion can be made from this spreadsheet, and it is given for general information only. For more information, please see your local bioinformatician and/or open the relevant Kraken report.");
  open(KRAKENREPORT,">","$outdir/report.tsv") or die "ERROR: could not open $outdir/report.tsv for writing: $!";
  print KRAKENREPORT join("\n",@report)."\n";
  close KRAKENREPORT;

  return $outdir;
}

sub runKraken{
  my($sample,$sampledir,$settings)=@_;

  my $html="$sampledir/report.html";
  return if(-e $html);

  my $reads=join(" ",@{ $$sample{fastq} });
  return if(!$reads);
  
  command("$KRAKENDIR/kraken --fastq-input --paired --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out $reads");

  command("$KRAKENDIR/kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  command("$KRAKENDIR/kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");

  command("$KRONADIR/ktImportText -o $html $sampledir/kraken.taxonomy");
}

sub reportContamination{
  my($sampledir,$expectedSpecies,$settings)=@_;
  return 100 if(!$expectedSpecies);

  my $taxfile="$sampledir/kraken.taxonomy";

  my %bestGuess;
  my $numCorrectReads=0;
  my $numContaminantReads=0;
  open(TAXONOMY,'<',$taxfile) or die "ERROR: could not open $taxfile for reading: $!";
  while(<TAXONOMY>){
    chomp;
    my($numReads,undef, $domain, $kingdom, $phylum, $class, $order, $family, $genus, $species)=split /\t/;
    $genus||="";
    $species||="";
    $species=~s/^.+\s+(.+)/$1/; # sometimes there are genus and species in the species column, but you just want the second word to be the species.
    my $scientificName="$genus $species";

    if($expectedSpecies eq $genus || $expectedSpecies eq $species || $expectedSpecies eq $scientificName){
      $numCorrectReads+=$numReads;
    } else {
      $numContaminantReads+=$numReads;
    }

    # Decide on a best guess for what this taxon is.
    for my $tier($scientificName, $genus, $family, $order, $class, $phylum, $kingdom, $domain){
      $tier||="";
      next if($tier=~/^\s*$/); # don't consider this tier if it's empty
      $bestGuess{$tier}+=$numReads;
      last;
    }
  }
  close TAXONOMY;

  my $bestGuess=(sort{$bestGuess{$b} <=> $bestGuess{$a}} keys(%bestGuess))[0];
  
  my $percentContamination=$numContaminantReads/($numContaminantReads+$numCorrectReads) * 100;
  return ($percentContamination, "$sampledir/report.html", $bestGuess) if wantarray;
  return $percentContamination;
}

sub usage{
  "Finds contamination in a miseq run
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  "
}

