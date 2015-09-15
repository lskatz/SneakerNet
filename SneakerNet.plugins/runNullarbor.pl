#!/usr/bin/env perl
# Runs nullarbor

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

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig logmsg samplesheetInfo command/;

$ENV{PATH}="$ENV{PATH}:/opt/mlst-1.2/bin:/opt/nullarbor/bin:/opt/abricate/bin:/opt/bcftools-1.2:/opt/samtools-1.2:/opt/megahit-1.0.2/bin:/opt/prokka-1.11/bin:/opt/snippy-2.6/bin:/opt/kraken";
$ENV{KRAKEN_DB_PATH}="/opt/kraken/minikraken_20141208";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug test numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  
  my $dir=$ARGV[0];

  nullarbor($dir,$settings);

  return 0;
}

sub nullarbor{
  my($dir,$settings)=@_;
  my $allsamples=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  # Run nullarbor by species, and so we need species names
  my %species;
  while(my($samplename,$info)=each(%$allsamples)){
    next if(ref($info) ne 'HASH');
    $species{$$info{species}}=1;
  }
  my @species=keys(%species);

  for my $species(@species){
    nullarborBySpecies($dir,$species,$settings);
  }
}

sub nullarborBySpecies{
  my($dir,$species,$settings)=@_;
  
  my $tsv=makeTabFile($dir,$species,$settings);
  my $mlstScheme=chooseMlstScheme($species,$settings);
  my $ref=chooseRef($dir,$species,$settings);

  my $outdir="$dir/$species.nullarbor";
     $outdir=~s/\s+/_/g;
  
  command("nullarbor.pl --name $species --mlst $mlstScheme --ref $ref --input $tsv --outdir $outdir --cpus $$settings{numcpus} 2>&1 | tee $outdir.log");
  command("nice make -C $outdir");

  return $outdir;
}

sub makeTabFile{
  my($dir,$species,$settings)=@_;
  $species||="all";

  my $allsamples=samplesheetInfo("$dir/SampleSheet.csv",$settings);
  my $tsv="$dir/samples.$species.tsv";
  open(TAB,">",$tsv) or die "ERROR: could not open $tsv for writing: $!";
  while(my($samplename,$info)=each(%$allsamples)){
    next if(ref($info) ne 'HASH');
    if($species eq 'all' || $$info{species} eq $species){
      print TAB join("\t",$samplename,@{$$info{fastq}})."\n";
    }
  }
  close TAB;
  return $tsv;
}

sub chooseMlstScheme{
  my($species,$settings)=@_;
  
  # TODO put this logic into a config file
  my $scheme="";
  if($species=~/listeria|monocytogenes/i){
    $scheme="lmonocytogenes";
  } elsif($species=~/Escherichia|E\.coli/i){
    $scheme="ecoli";
  } elsif($species=~/Campy/i){
    $scheme="campylobacter";
  } elsif($species=~/cholerae/i){
    $scheme="vcholerae";
  } elsif($species=~/vibrio/i){
    $scheme="vibrio";
  } elsif($species=~/salmonella/i){
    $scheme="senterica";
  } else {
    $scheme="";
  }
  return $scheme;
}

sub chooseRef{
  my($dir,$species,$settings)=@_;
  my $ref="$dir/ref.fasta";
  open(REF,">",$ref) or die "ERROR: could not open $ref for writing: $!";
  print REF ">BLANKREF\nAAAAAAAAAAAAAAAAAAAAAAAAA\nTTTTTTTTTTTTTTTTTTT\n";
  close REF;
  return $ref;
}
################
# Utility subs #
################

sub usage{
  "Runs nullarbor for the read set
  Usage: $0 runDir
  "
}

