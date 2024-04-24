#!/usr/bin/env perl
# Use Kraken to find contamination.
# Modeling the workflow from https://github.com/lskatz/lskScripts/blob/205c03d3c4c44bfa2f961015c8d5f6971e5698d2/qsub/launch_kraken.sh

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/cp mv/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "5.0";
our $CITATION= "Detect contamination with Kraken plugin by Lee Katz.  Uses Kraken1.";

# Get the executable directories
my $tmpSettings=readConfig();
my @contaminatedSamples; # list of potentially contaminated samples
my $warningsMsg = ""; # for warnings on the html report
my $errorsMsg   = ""; # for errors on the html report

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(minpercent|min_percent|min-percent=f minwarning|min_warning|min-warning=f version citation check-dependencies help debug tempdir=s numcpus=i force)) or die $!;
  my @exe = qw(zip kraken kraken-translate kraken-report ktImportText);
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe       => \@exe,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{KRAKEN_DEFAULT_DB} ||= die "ERROR: KRAKEN_DEFAULT_DB needs to be defined under config/settings.conf";
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);
  $$settings{minpercent} ||= 25.0;
  $$settings{minwarning} ||= 5.0;

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";
  mkdir "$dir/SneakerNet/kraken";
  
  my $outdir=runKrakenOnDir($dir,$settings);

  # make the report emailable 
  cp("$outdir/report.tsv", "$dir/SneakerNet/forEmail/kraken.tsv");
  #command("cd $dir/SneakerNet/forEmail && zip -v9 kraken.zip *.kraken.html && rm -v *.kraken.html");

  if(@contaminatedSamples){
    $warningsMsg .= "There is at least one potentially contaminated sample ( >=$$settings{minwarning}%) ";
  }

  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{
    version        => $VERSION,
    krakenDatabase => $$settings{KRAKEN_DEFAULT_DB},
    table          => "$dir/SneakerNet/forEmail/kraken.tsv",
    warnings       => $warningsMsg,
    mqc            => $rawMultiQC,
    exe            => \@exe,
  });
  # also record if there are any potentially contaminated samples

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/kraken.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/kraken_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Contamination detection (Kraken)\"\n";
  print $outFh "#description: \"Estimating contamination using Kraken1.<br />Database $$settings{KRAKEN_DEFAULT_DB}<br />$plugin v$VERSION $docLink $pluginLink\"\n";
  print $outFh "#anchor: '$anchor'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  while(<$fh>){
    next if(/^#/);
    print $outFh $_;
  }
  close $fh;

  return $outtable;
}

sub runKrakenOnDir{
  my($dir,$settings)=@_;
  my $outdir="$dir/SneakerNet/kraken";
  system("mkdir -p $outdir");

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my %filesToTransfer=(); # hash keys are species names
  my @report; # reporting contamination in an array, in case I want to sort it later
  push(@report, join("\t", qw(
    NAME LABELED_TAXON 
    BEST_GUESS PERCENTAGE_OF_GENOME_IS_BEST_GUESS 
    RANK
    MAJOR_CONTAMINANT PERCENTAGE_CONTAMINANT
  )));
  while(my($sampleName,$s)=each(%$sampleInfo)){

    my $sampledir="$outdir/$sampleName";
    system("mkdir -p $sampledir");
    logmsg "Running Kraken on $sampleName";
    logmsg "  Database: $$settings{KRAKEN_DEFAULT_DB}";

    my $expectedSpecies=$$s{species} || $$s{taxon} || "UNKNOWN";
    #my($percentTaxon,$html,$bestGuess)=guessTaxon($sampledir,$settings);

    # Create hash of bestGuess
    # e.g., {guess=>\%guessedTaxon, contaminant=>\%majorConflictingTaxon, html=>"$sampledir/report.html"};
    my $guesses = guessTaxon($sampledir,$settings);
    my $html = $$guesses{html};
    if(!$html){
      logmsg "WARNING: html file not defined for $sampleName";
      next;
    }
    if(!-e $html){
      logmsg "WARNING: html file not found: $html";
      next;
    }

    # Raise a warning if there is a significant amount of contamination
    if($$guesses{contaminant}{percent} >= $$settings{minwarning}){
      raiseWarning($sampleName, $settings);
    }

    # Add onto the contamination report
    push(@report,join("\t",
        $sampleName,$expectedSpecies,
        $$guesses{guess}{taxname}, $$guesses{guess}{percent},
        $$guesses{guess}{rank},
        $$guesses{contaminant}{taxname}, $$guesses{contaminant}{percent},
    ));
    logmsg "Including for email: $html";
    #cp($html,"$dir/SneakerNet/forEmail/$sampleName.kraken.html");

    # Report anything with >10% contamination to the printout.
    #if($percentContaminated > 10){
      #logmsg "$sampleName (taxon: $expectedSpecies) is $percentContaminated% contaminated";
    #}
  }

  # print the report to a file.
  # Note: Taylor reports 20% unclassified and a 4% contamination in one genome, which was considered uncontaminated.  Further testing is needed.
  push(@report,"# A genome with >90% reads in agreement with its species has not been shown to indicate contamination");
  push(@report,"# More testing needs to be performed before any conclusion can be made from this spreadsheet, and it is given for general information only. For more information, please see your local bioinformatician and/or open the relevant Kraken report.");
  open(KRAKENREPORT,">","$outdir/report.tsv") or die "ERROR: could not open $outdir/report.tsv for writing: $!";
  print KRAKENREPORT join("\n",@report)."\n";
  close KRAKENREPORT;

  return $outdir;
}

sub guessTaxon{
  my($sampledir,$settings)=@_;
  logmsg $sampledir;

  my $taxfile="$sampledir/kraken.report";

  # If sn_kraken did not complete, a file will not be present
  if(!-e $taxfile){
    logmsg "WARNING: kraken report not found at $taxfile. It is possible that sn_kraken.pl or kraken was not run on this sample yet.";
    return {}; # return empty hash because that is the var type expected
  }

  my %bestGuess;
  my @header = qw(percent classifiedReads specificClassifiedReads rank taxid taxname);
  open(TAXONOMY, '<', $taxfile) or die "ERROR: could not open $taxfile for reading: $!";
  while(my $taxline = <TAXONOMY>){
    # trim whitespace
    $taxline =~ s/^\s+|\s+$//g;

    # Split fields on tab, but also remove whitespace on fields
    my @F = split(/\s*\t\s*/, $taxline);

    # Name the fields
    my %field;
    @field{@header} = @F;

    # We only care about named ranks like phylum, family, genus, or species
    next if($field{rank} eq '-');

    push(@{ $bestGuess{$field{rank}} }, \%field);

  }
  close TAXONOMY;

  my @sortedRank = qw(S G F O C P K D U);

  # Starting with species, if >min-percent% reads are attributed to
  # a taxon, then guess that taxon.
  my %guessedTaxon;
  RANK:
  for my $rank(reverse @sortedRank){
    $bestGuess{$rank} //= [];
    for my $taxHash(sort {$$b{percent} <=> $$a{percent}} @{ $bestGuess{$rank} }){
      if($$taxHash{percent} > $$settings{minpercent}){
        %guessedTaxon = %$taxHash;
      }
    }
  }
  my $rank = $guessedTaxon{rank};

  # Are there any competing taxa at that rank?
  # Initialize the competing taxon to boolean false numbers
  my %majorConflictingTaxon = (rank=>"", percent => 0, taxid=>0, taxname=>".", specificClassifiedReads=>0, classifiedReads=>0);
  $bestGuess{$rank} //= [];
  for my $taxHash(sort {$$b{percent} <=> $$a{percent}} @{ $bestGuess{$rank} }){
    # Skip comparing against self, by comparing taxids
    next if($$taxHash{taxid} == $guessedTaxon{taxid});

    # Accept the conflicting taxon if it's over 1% and if
    # it is the biggest contaminant.
    if($$taxHash{percent} > 1 && $$taxHash{percent} > $majorConflictingTaxon{percent}){
      %majorConflictingTaxon = %$taxHash;
    }
  }

  # Return the best guess and the major conflict
  return {guess=>\%guessedTaxon, contaminant=>\%majorConflictingTaxon, html=>"$sampledir/report.html"};
}

sub raiseWarning{
  my($sampleName, $settings) = @_;
  push(@contaminatedSamples, $sampleName);
}

sub usage{
  print "Finds contamination in a miseq run with kraken and guesses the taxon
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
  --force         Overwrite any previous results
  --min-percent   What percent of reads have to be attributed
                  to a taxon before presuming it as the taxon
                  sequenced? Default: 25.0
  --min-warning   What percent do reads of a contaminant taxon
                  have to be before a warning is triggered
                  in the report? Default: 5
";
  exit(0);
}

