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

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.0";
our $CITATION = "Save failed genomes by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation help force tempdir=s debug numcpus=i coverage=f)) or die $!;
  $$settings{coverage} //= 10;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }
  if($$settings{citation}){
    print $CITATION."\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

  if(! -d "$dir/SneakerNet"){
    mkdir "$dir/SneakerNet";
  }
  if(! -d "$dir/SneakerNet/forEmail"){
    mkdir "$dir/SneakerNet/forEmail";
  }

  my $table = saveGenomes($dir, $settings);

  recordProperties($dir,{version=>$VERSION, table=>$table});

  return 0;
}

sub saveGenomes{
  my($dir, $settings)=@_;

  # Get sample information
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  
  # Read the readMetrics file into %readMetrics
  open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
  my @header=split(/\t/,<READMETRICS>); chomp(@header);
  my %readMetrics = ();
  while(<READMETRICS>){
    chomp;
    my %F;
    @F{@header}=split(/\t/,$_);
    $F{File} = basename($F{File});
    $readMetrics{$F{File}} = \%F;
  }
  close READMETRICS;

  # Calculate coverage for each sample
  my %saved = ();
  for my $samplename(keys(%$sampleInfo)){
    my $totalCoverage = 0;
    
    # Get the name of all files linked to this file through the sample.
    $$sampleInfo{$samplename}{fastq}//=[]; # Set {fastq} to an empty list if it does not exist
    my @file=@{$$sampleInfo{$samplename}{fastq}};
    for my $fastq(@file){
      my $fastqMetrics = $readMetrics{basename($fastq)};

      # Coverage
      if($$fastqMetrics{coverage} eq '.'){ # dot means coverage is unknown
        $totalCoverage = -1; # -1 means 'unknown' coverage
        logmsg "$fastq unknown coverage" if($$settings{debug});
      } else {
        $$fastqMetrics{coverage} ||= 0; # force it to be a number if it isn't already
        $totalCoverage += $$fastqMetrics{coverage};
        logmsg "Sample $samplename += $$fastqMetrics{coverage}x => ${totalCoverage}x" if($$settings{debug});
      }
    }

    if(
         $totalCoverage >= $$settings{coverage} 
      && $totalCoverage <  $$sampleInfo{$samplename}{taxonRules}{coverage}
    ){
      $saved{$samplename} = \@file;
    }
  }

  # Save the samples with an rsync
  for my $samplename(keys(%saved)){
    rsync($$sampleInfo{$samplename}, $settings);
  }

  # Make a table of results
  my $table = "sample\tfastqs\n";
  while(my($name, $fastqs) = each(%saved)){
    $table .= $name;
    $table .= "\t".join(";", @$fastqs);
    $table .= "\n";
  }
  my $outfile = "$dir/SneakerNet/forEmail/qc_fails.tsv";
  open(my $fh, '>', $outfile) or die "ERROR writing to $outfile: $!";
  print $fh $table;
  close $fh;

  return $outfile;
}

sub rsync{
  my($sample, $settings) = @_;

  my $subfolder = "$$sample{taxonRules}{dest_subfolder}/QC_Fails";
  my $fileString = join(" ", @{$$sample{fastq}});

  command("rsync -av -q --no-motd -av --no-g --copy-links --chmod=Du=rwx,Dg=rx,Do=rx,Fu=rw,Fg=r,Fo=r $fileString $$settings{transfer_destination_string}/$subfolder/");

  return 1;
}
    

sub usage{
  print "Saves genomes into \$subfolder/QC_Fails if they
  do not meet coverage but meet a lower coverage threshold.

  Usage: $0 MiSeq_run_dir
  --coverage  10  Minimum coverage to save a genome, but
                  the coverage found in
                  conf/taxonProperties.conf (usually at
                  least 20x)
  \n";
  exit 0;
}

