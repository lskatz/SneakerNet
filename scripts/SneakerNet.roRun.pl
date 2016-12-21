#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i outdir=s)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  $$settings{numcpus}||=1;
  $$settings{outdir}||="sneakernet.out";

  die usage() if($$settings{help} || !@ARGV);

  for my $dir(@ARGV){
    my $sneakernetDir = makeSneakernetDir($dir,$settings);
    saveSneakernetDir($sneakernetDir, $$settings{outdir});
  }

  return 0;
}

sub makeSneakernetDir{
  my($dir,$settings)=@_;

  my $outdir="$$settings{tempdir}/runData";
  mkdir $outdir;

  my @fastq       = glob("$dir/Data/Intensities/BaseCalls/*.fastq.gz");
  my $sampleSheet=  "$dir/Data/Intensities/BaseCalls/SampleSheet.csv";
  my $config      =  "$dir/Data/Intensities/BaseCalls/config.xml";
  my @interop     =  glob("$dir/InterOp/*");
  my @xml         = ("$dir/CompletedJobInfo.xml",
                     "$dir/runParameters.xml",
                     "$dir/GenerateFASTQRunStatistics.xml",
                     "$dir/RunInfo.xml",
                    );
  
  if(!@fastq){
    logmsg "WARNING: no fastq files were found; attempting bcl2fastq";
    @fastq=bcl2fastq($dir,$settings);
    if(!@fastq){
      die "ERROR: could not find any fastq files in $dir";
    }
  }
  if(!@interop){
    die "ERROR: no interop files were found in $dir";
  }
  for(@fastq,$sampleSheet, $config, @interop, @xml){
    if(!-e $_){
      die "ERROR: file does not exist: $_";
    }
  }

  for(@fastq, $sampleSheet, $config){
    my $to="$outdir/".basename($_);
    cp($_, $to);
  }
  mkdir "$outdir/QC";
  for(@xml){
    my $to="$outdir/QC/".basename($_);
    cp($_, $to);
  }
  mkdir "$outdir/QC/InterOp";
  for(@interop){
    my $to="$outdir/QC/InterOp/".basename($_);
    cp($_, $to);
  }

  return $outdir;
}

sub bcl2fastq{
  my($dir,$settings)=@_;
  my $fastqdir="$$settings{tempdir}/bcl2fastq";
  mkdir($fastqdir);

  command("bcl2fastq --input-dir $dir/Data/Intensities/BaseCalls --runfolder-dir $dir --output-dir $fastqdir --processing-threads $$settings{numcpus} --demultiplexing-threads $$settings{numcpus} --barcode-mismatches 1 >&2");

  my @fastq=glob("$$settings{tempdir}/bcl2fastq/*.fastq.gz");
  return @fastq;
}

sub saveSneakernetDir{
  my($tmpdir,$outdir,$settings)=@_;
  File::Copy::mv($tmpdir,$outdir) or die "ERROR: could not move $tmpdir to $outdir: $!";
  return 1;
}

sub cp{
  my($from,$to)=@_;
  if(-e $to && -s $to > 0){
    logmsg "Found $to. Not copying";
    return 1;
  }
  logmsg "cp $from to $to";
  my $return=File::Copy::cp($from,$to) or die "ERROR: could not copy $from to $to: $!";
  return $return;
}

sub usage{
  "Parses an unaltered Illumina run and formats it
  into something usable for SneakerNet

  Usage: $0 illuminaDirectory [illuminaDirectory2...]
  
  --numcpus  1
  "
}
