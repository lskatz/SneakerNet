#!/usr/bin/env perl
# Run wgMLST (or cgMLST or whatever genome-wide MLST scheme)

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/cp mv/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use Cwd qw/realpath/;
use FindBin qw/$RealBin/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "1.2";
our $CITATION= "wgMLST plugin by Lee Katz.  Uses chewBBACA for wgMLST.";

# wget --recursive http://enterobase.warwick.ac.uk/schemes/SALwgMLST.cgMLSTv1/
#                  http://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/
local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help debug tempdir=s numcpus=i force )) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      zip       => 'zip --version | grep "This is Zip"',
      'python3 (Python3)'   => 'python3 --version',
      'chewBBACA.py (chewBBACA)' => 'chewBBACA.py -h 2>&1 | grep -m 1 -i version',
      'blastn (BLAST+)'    => 'blastn -version | head -n 1',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];

  my $outdir = "$dir/SneakerNet/wgMLST";
  mkdir $outdir;

  runWgMlst($dir, $settings);

  # Let's make an output file
  my $outputTsv = "$dir/SneakerNet/forEmail/wgMLST.tsv";
  my @header = qw(Genome db  EXC     INF     LNF     PLOT    NIPH    ALM     ASM);
  open(my $fh, '>', $outputTsv) or die "ERROR: could not write to $outputTsv: $!";
  print $fh join("\t", @header)."\n";
  for my $sampleDir(glob("$outdir/*")){
    next if(!-d $sampleDir);
    my $summaryIn = "$sampleDir/results_statistics.tsv";
    open(my $inFh, '<', $summaryIn) or die "ERROR: could not read from $summaryIn: $!";

    # Get this sample's header
    my $summaryHeader = <$inFh>;
    chomp($summaryHeader);
    my @summaryHeader = split(/\t/, $summaryHeader);

    # Get the result for this sample
    my $res = <$inFh>;
    chomp($res);
    close $inFh;

    # get each value attached to its proper header
    my %F;
    my @F = split(/\t/, $res);
    @F{@summaryHeader} = @F;

    # Additional information about this run
    open(my $optsFh, '<', "$sampleDir/sneakernet.tsv") or die "ERROR: could not read from $sampleDir/sneakernet.tsv: $!";
    while(<$optsFh>){
      chomp;
      my($key,$value) = split(/\t/);
      $F{lc($key)}=$value;
    }

    for my $h(@header){
      print $fh "$F{$h}\t";
    }
    print $fh "\n";
  }
  print $fh "# EXC: alleles which have exact matches (100% DNA identity) with previously identified alleles\n";
  print $fh "# INF: inferred new alleles using Prodigal CDS predictions\n";
  print $fh "# LNF: loci not found.\n";
  print $fh "# PLOT: possible loci on the tip of the query genome contigs.\n";
  print $fh "# NIPH: non-informative paralogous hit\n";
  print $fh "# NIPHEM: similar to NIPH classification (NIPH with exact match), but specifically referring to exact matches\n";
  print $fh "# ALM: alleles 20% larger than length mode of the distribution of the matched loci\n";
  print $fh "# ASM: similar to ALM but for alleles 20% smaller than length mode distribution of the matched loci\n";
  close $fh;

  recordProperties($dir,{version=>$VERSION,table=>$outputTsv});

  return 0;
}

sub runWgMlst{
  my($dir,$settings)=@_;

  my $outdir="$dir/SneakerNet/wgMLST";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    my $sampleOutdir  = "$dir/SneakerNet/wgMLST/$sampleName";
    if(-d $sampleOutdir && !$$settings{force}){
      logmsg "SKIP: found $sampleOutdir";
      next;
    }
    # If the sample output directory is still there (ie,
    # --force), then remove the target first.
    if(-d $sampleOutdir){
      command("rm -rf '$sampleOutdir'");
    }
    
    # Grab the assembly
    my @asm = glob("$dir/SneakerNet/assemblies/$sampleName/*.fasta");
    my $asm = $asm[0];
    if(!$asm || !-e $asm){
      logmsg "WARNING: no assembly found for $sampleName. Not analyzing";
      next;
    }
    # Assembly needs to be > 20000bp to predict genes from.
    # Be safe with 30000 bytes.
    if(-s $asm < 30000){
      logmsg "WARNING: assembly is too small for $sampleName. Not analyzing.";
      next;
    }

    # Database directory, formatted by chewBBACA
    my $wgMlstSubdir = $$s{taxonRules}{wgMLST};
    if(!defined($wgMlstSubdir)){
      logmsg "No wgMLST scheme defined in taxonProperties.conf for $sampleName. Skipping";
      logmsg "Taxon was defined as ".$$s{taxon};
      next;
    }

    logmsg "Taxon was defined for $sampleName. prepping the analysis.";
    my $wgMlstDir    = realpath($RealBin."/../db/wgMLST/$wgMlstSubdir");
    my $wgMlstDirBak = $wgMlstDir;
    # Error checking for the wgMLST database: does the dir exist
    if(!-d $wgMlstDir){
      logmsg "WARNING: no wgMLST directory found in $wgMlstDir - not analyzing";
      logmsg "To download a scheme, follow the instructions at https://github.com/Public-Health-Bioinformatics/pubmlst_client";
      next;
    }
    my @allLoci      = glob("$wgMlstDir/*.fasta");
    my @allShortLoci = glob("$wgMlstDir/short/*.fasta");
    my $numLoci      = @allLoci;
    my $numShortLoci = @allShortLoci;
    # Error checking for the wgMLST database: does the dir have loci
    if(!-d "$wgMlstDir/short" || $numLoci==0 || $numShortLoci==0){
      logmsg "WARNING: the database might not be formatted under $wgMlstDir - please see chewBBACA documentation on how to format the database properly.";
      next;
    }
    # Error checking for the wgMLST database: do the number of loci
    # in the short directory equal to the number of loci overall
    if($numLoci != $numShortLoci){
      logmsg "WARNING: the database at $wgMlstDir might not have finished formatting but I will analyze anyway";
    }

    # Results to a temporary directory
    my $tmpout = "$$settings{tempdir}/mlst-wg-$sampleName";
    mkdir $tmpout;

    # Run chewBBACA
    # Tests show 3.5 minutes on 12 cpus for 224 loci
    if(-d "/dev/shm"){
      # Copy over the database to ram since it is disk IO intensive
      my $ramTempdir = "/dev/shm/$ENV{USER}";
      mkdir $ramTempdir;
      my $newWgMlstDir = tempdir("$0.$$s{taxon}.XXXXXX",DIR=>$ramTempdir, CLEANUP=>1);
      logmsg "Copying wgMLST database to RAM at $newWgMlstDir";
      system("cp -r $wgMlstDir/* $newWgMlstDir/ >&2");
      $wgMlstDir = $newWgMlstDir;
    }
    logmsg "Running chewBBACA.py AlleleCall on $asm on scheme $wgMlstDir ($numLoci loci)";
    command("chewBBACA.py AlleleCall --fr -i $asm -g $wgMlstDir --cpu $$settings{numcpus} -o $tmpout");

    # move over the results
    logmsg "wgMLST results in $tmpout";
    my @chewbbacaOut = glob("$tmpout/results_*");
    my $from = $chewbbacaOut[0];
    if(!defined($from)){
      logmsg "ERROR: no results found!";
      logmsg `tree $tmpout`;
      next;
    }
    logmsg "Found results in $from";
    command("cp -r '$from' '$sampleOutdir'");

    my $sneakerNetSummary = "$sampleOutdir/sneakernet.tsv";
    open(my $fh, '>', $sneakerNetSummary) or die "ERROR: could not write to $sneakerNetSummary: $!";
    print $fh "db\t$wgMlstDirBak\n";
    close $fh;
  }

  return "$dir/SneakerNet/wgMLST";
}

sub usage{
  print "Runs wgMLST
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --force
";
  exit(0);
}

