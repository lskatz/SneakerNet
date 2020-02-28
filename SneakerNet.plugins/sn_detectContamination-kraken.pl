#!/usr/bin/env perl
# Use Kraken to find contamination.
# Modeling the workflow from https://github.com/lskatz/lskScripts/blob/205c03d3c4c44bfa2f961015c8d5f6971e5698d2/qsub/launch_kraken.sh

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/cp/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "2.2";
our $CITATION= "Detect contamination with Kraken plugin by Lee Katz.  Uses Kraken1.";

# Get the executable directories
my $tmpSettings=readConfig();
#my $KRAKENDIR=$$tmpSettings{KRAKENDIR} || die "ERROR: could not find KRAKENDIR in config";
#my $KRONADIR=$$tmpSettings{KRONADIR} || die "ERROR: could not find KRONADIR in config";
#$ENV{PATH}="$ENV{PATH}:$KRAKENDIR:$KRONADIR";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(minpercent|min_percent|min-percent=f version citation check-dependencies help debug tempdir=s numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      zip       => 'zip --version | grep "This is Zip"',
      kraken    => 'kraken --version | grep -m 1 version',
      'kraken-translate' => 'kraken-translate --version | grep -m 1 version',
      'kraken-report'    => 'kraken-report --version | grep -m 1 version',
      'ktImportText'     => 'ktImportText | grep "/" | grep -P -m 1 -o "KronaTools .*ktImportText"',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{KRAKEN_DEFAULT_DB} ||= die "ERROR: KRAKEN_DEFAULT_DB needs to be defined under config/settings.conf";
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);
  $$settings{minpercent} ||= 25;

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";
  mkdir "$dir/SneakerNet/kraken";
  
  my $outdir=runKrakenOnDir($dir,$settings);

  # make the report emailable 
  cp("$outdir/report.tsv", "$dir/SneakerNet/forEmail/kraken.tsv");
  command("cd $dir/SneakerNet/forEmail && zip -v9 kraken.zip *.kraken.html && rm -v *.kraken.html");

  recordProperties($dir,{version=>$VERSION,krakenDatabase=>$$settings{KRAKEN_DEFAULT_DB},table=>"$dir/SneakerNet/forEmail/kraken.tsv"});

  return 0;
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
    logmsg "  Database: $$settings{KRAKEN_DEFAULT_DB}";
    my $krakenWorked=runKraken($s,$sampledir,$settings);

    if(!$krakenWorked){
      logmsg "Kraken was not completed successfully on $sampleName. I will not output results for $sampleName";
      next;
    }

    my $expectedSpecies=$$s{species} || $$s{taxon} || "UNKNOWN";
    #my($percentTaxon,$html,$bestGuess)=guessTaxon($sampledir,$settings);

    # Create hash of bestGuess
    # e.g., {guess=>\%guessedTaxon, contaminant=>\%majorConflictingTaxon, html=>"$sampledir/report.html"};
    my $guesses = guessTaxon($sampledir,$settings);
    my $html = $$guesses{html};
    if(!$html){
      logmsg "WARNING: html file not found: $html";
      next;
    }

    # Add onto the contamination report
    push(@report,join("\t",
        $sampleName,$expectedSpecies,
        $$guesses{guess}{taxname}, $$guesses{guess}{percent},
        $$guesses{guess}{rank},
        $$guesses{contaminant}{taxname}, $$guesses{contaminant}{percent},
    ));
    logmsg "Including for email: $html";
    cp($html,"$dir/SneakerNet/forEmail/$sampleName.kraken.html");

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

sub runKraken{
  my($sample,$sampledir,$settings)=@_;

  my $html="$sampledir/report.html";
  if(-e $html && !$$settings{force}){
    return 1;
  }

  if(!defined($$sample{fastq})){
    logmsg "ERROR: no reads found for $sampledir";
    return 0;
  }
  
  # Skip small file sizes.
  # TODO: use something better like readMetrics.pl 
  for(@{ $$sample{fastq} }){
    if(-s $_ < 10000){
      logmsg "There are few reads in $$sample{sample_id}. Skipping.";
      return 0;
    }
  }
  
  # Force an array
  if(ref($$sample{fastq}) ne 'ARRAY'){
    $$sample{fastq} = [$$sample{fastq}];
  }
  my @fastq = @{$$sample{fastq}};

  if(@fastq == 1){
    return runKrakenSE(@_);
  }
  elsif(@fastq == 2){
    return runKrakenPE(@_);
  } else {
    my %sampleCopy = %$sample;
    splice(@{ $sampleCopy{fastq} }, 2);
    return runKrakenPE(\%sampleCopy, $sampledir, $settings);
  }

  logmsg "INTERNAL ERROR";
  return 0;
}
sub runKrakenSE{
  my($sample,$sampledir,$settings)=@_;
  my $html="$sampledir/report.html";

  my $reads = $$sample{fastq}[0];

  return 0 if(!$reads);

  command("kraken --fastq-input $reads --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out ");

  command("kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  command("kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");

  # To capture unclassified reads, we can get the third
  # column of the first row of the report file. This
  # information can be appended to the taxonomy file
  # on the last line.
  open(my $reportFh, "<", "$sampledir/kraken.report") or die "ERROR: could not read $sampledir/kraken.report: $!";
  my $firstLine=<$reportFh>;
  close $reportFh;
  my $unclassifiedReadsCount=(split(/\t/, $firstLine))[2];
  open(my $taxFh, ">>", "$sampledir/kraken.taxonomy") or die "ERROR: could not append to $sampledir/kraken.taxonomy: $!";
  print $taxFh $unclassifiedReadsCount."\n";
  close $taxFh;

  command("ktImportText -o $html $sampledir/kraken.taxonomy");

  # Go ahead and remove kraken.out which is a huge file
  unlink("$sampledir/kraken.out");

  if(! -e "$sampledir/kraken.taxonomy"){
    return 0;
  }

  return 1;
}

sub runKrakenPE{
  my($sample,$sampledir,$settings)=@_;
  my $html="$sampledir/report.html";

  my @twoReads = (@{$$sample{fastq}})[0,1];
  my $reads="'".join("' '", @twoReads)."'";
  return 0 if(!$reads);
  
  command("kraken --fastq-input --paired $reads --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out ");

  command("kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  command("kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");

  # To capture unclassified reads, we can get the third
  # column of the first row of the report file. This
  # information can be appended to the taxonomy file
  # on the last line.
  open(my $reportFh, "<", "$sampledir/kraken.report") or die "ERROR: could not read $sampledir/kraken.report: $!";
  my $firstLine=<$reportFh>;
  close $reportFh;
  my $unclassifiedReadsCount=(split(/\t/, $firstLine))[2];
  open(my $taxFh, ">>", "$sampledir/kraken.taxonomy") or die "ERROR: could not append to $sampledir/kraken.taxonomy: $!";
  print $taxFh $unclassifiedReadsCount."\n";
  close $taxFh;

  command("ktImportText -o $html $sampledir/kraken.taxonomy");

  # Go ahead and remove kraken.out which is a huge file
  unlink("$sampledir/kraken.out");

  if(! -e "$sampledir/kraken.taxonomy"){
    return 0;
  }

  return 1;
}

sub guessTaxon{
  my($sampledir,$settings)=@_;
  logmsg $sampledir;

  my $taxfile="$sampledir/kraken.report";

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
  for my $rank(@sortedRank){
    for my $taxHash(sort {$$b{percent} <=> $$a{percent}} @{ $bestGuess{$rank} }){
      if($$taxHash{percent} > $$settings{minpercent}){
        %guessedTaxon = %$taxHash;
        last RANK;
      }
    }
  }
  my $rank = $guessedTaxon{rank};

  # Are there any competing taxa at that rank?
  # Initialize the competing taxon to boolean false numbers
  my %majorConflictingTaxon = (rank=>"", percent => 0, taxid=>0, taxname=>".", specificClassifiedReads=>0, classifiedReads=>0);
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

sub usage{
  print "Finds contamination in a miseq run with kraken and guesses the taxon
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
  --force         Overwrite any previous results
  --min-percent   What percent of reads have to be attributed
                  to a taxon before presuming it as the taxon
                  sequenced? Default: 25
";
  exit(0);
}

