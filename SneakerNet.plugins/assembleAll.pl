#!/usr/bin/env perl
# Transfers files to a remote destination and QCs them beforehand.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use FindBin;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig samplesheetInfo command logmsg fullPathToExec/;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help tempdir=s debug numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  # Check for required executables
  for (qw(megahit run_assembly_metrics.pl)){
    fullPathToExec($_);
  }
 
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/assemblies";
  assembleAll($dir,$settings);

  return 0;
}

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");

    my $assembly=assembleSample($sample,$info,$settings);
    next if(!$assembly);
    my $outdir="$dir/SneakerNet/assemblies/$sample";
    mkdir $outdir;
    cp($assembly,$outdir) or die "ERROR copying $assembly to $outdir/: $!";

    # TODO run prokka
  }
  
  # TODO run assembly metrics with min contig size=0.5kb
  # TODO compare genome size vs expected
}

sub assembleSample{
  my($sample,$sampleInfo,$settings)=@_;

  my $R1=$$sampleInfo{fastq}[0];
  my $R2=$$sampleInfo{fastq}[1];
  if(!$R1){
    logmsg "Could not find R1 for $sample. Skipping";
    return "";
  }
  if(!$R2){
    logmsg "Could not find R2 for $sample. Skipping";
    return "";
  }

  logmsg "Assembling $sample";

  my $outdir="$$settings{tempdir}/$sample";
  
  system("rm -rf '$outdir'"); # make sure any previous runs are gone
  command("megahit -1 $R1 -2 $R2 --out-dir '$outdir' -t $$settings{numcpus} 1>&2");
  die "ERROR with running megahit on $sample: $!" if $?;

  return "$outdir/final.contigs.fa";
}

sub usage{
  "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  "
}

