#!/usr/bin/env perl
# Hello World example

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.0";
our $CITATION = "Genotyping SneakerNet plugin by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      kma       => "kma -v",
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

  if(! -d "$dir/SneakerNet"){
    mkdir "$dir/SneakerNet";
  }
  if(! -d "$dir/SneakerNet/forEmail"){
    mkdir "$dir/SneakerNet/forEmail";
  }
  if(! -d "$dir/SneakerNet/genotype"){
    mkdir "$dir/SneakerNet/genotype";
  }

  my $table = genotypeSamples($dir, $settings);

  my $kmaVersion = `kma -v`;
  chomp($kmaVersion);
  recordProperties($dir,{version=>$VERSION, table=>$table, kma_version=>$kmaVersion});

  return 0;
}

sub genotypeSamples{
  my($dir, $settings) = @_;

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");

    # See if there is a genotype database
    my $genotypeDatabases = [];
    if($$info{taxonRules}{genotype}){
      # If there is a database for this taxon, then
      # check to make sure it has been downloaded.
      $genotypeDatabases = checkGenotypeDatabases($$info{taxonRules}{genotype}, $settings);
    }
    else {
      logmsg "No databases listed for sample $sample under taxon $$info{taxon}";
      logmsg "Therefore, I will not genotype $sample.";
      next;
    }

    # Genotype with the current databases
    my $outdir = "$dir/SneakerNet/genotype/$sample";
    mkdir($outdir);
    # Genotype one database at a time
    my $numDatabases = @$genotypeDatabases;
    for(my $i=0; $i<$numDatabases; $i++){
      my $R1 = $$info{fastq}[0];
      my $R2 = $$info{fastq}[1];
      my $target = "$outdir/".basename($$genotypeDatabases[$i]);
      if(! -e "$target.res" || $$settings{force}){
        command("kma -ipe $R1 $R2 -o $target -t_db $$genotypeDatabases[$i] -ID 60.0");
      }
    }
  }

  # Make an output summary table
  my $table = "$dir/SneakerNet/forEmail/genotype.tsv";
  command("head -n 1 -q $dir/SneakerNet/genotype/*/*.res | head -n 1 > $table");
  command("tail -n +2 -q $dir/SneakerNet/genotype/*/*.res | sort | uniq >> $table");
  return $table;
}
    
# Download any databases
sub checkGenotypeDatabases{
  my($urls, $settings) = @_;
  my @db;

  if(!-e "$RealBin/../db/genotype"){
    mkdir("$RealBin/../db/genotype");
  }

  $urls = (ref($urls) eq 'ARRAY')?$urls:[$urls];

  for my $u(@$urls){
    my $target = "$RealBin/../db/genotype/".basename($u);

    # Download if the fasta is not there yet
    if(!-e $target){
      command("wget '$u' -O $target");
    }

    # Index with kma
    if(!-e "$target.length.b"){
      command("kma index -i $target");
    }

    push(@db, $target);
  }

  return \@db;
}

sub usage{
  print "Run genotyping. Fasta databases can be listed in taxonProperties.conf.

  Usage: $0 MiSeq_run_dir
  --numcpus  1  Number of CPUs
  --force       Whether to overwrite results
  --debug       Print debugging information
  --version     Print the version and exit
  --help        This usage menu
  --tempdir  '' Specify the temporary directory
  ";
  exit 0;
}

