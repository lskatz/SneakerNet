#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig samplesheetInfo_tsv command logmsg version/;

our $VERSION = "1.0";

my @fastqExt=qw(.fastq.gz .fq.gz .fastq .fq);

my $snVersion=version();

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version help numcpus=i debug tempdir=s force)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

  mkdir("$dir/SneakerNet/SalmID");
  mkdir "$dir/SneakerNet/forEmail";
  
  identifyEach($dir,$settings);

  command("cat $dir/SneakerNet/SalmID/*/*.tsv > $dir/SneakerNet/forEmail/SalmID.tsv");

  return 0;
}

sub identifyEach{
  my($dir,$settings)=@_;

  my $sampleInfo=samplesheetInfo("$dir/samples.tsv",$settings);
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases
    my $taxon=$$s{species} || $$s{taxon} || 'NOT LISTED';

    my $outdir="$dir/SneakerNet/SalmID/$sampleName";
    next if(-e $outdir);
    my $tempdir="$$settings{tempdir}/$sampleName";
    mkdir $tempdir;
    for my $fastq(@{ $$s{fastq} }){
      logmsg $fastq;
      command("SalmID.py -i '$fastq' > $tempdir/".basename($fastq,@fastqExt).".tsv");
    }
    command("mv $tempdir $outdir");
  }
}

sub usage{
  "Run SalmID on each fastq file to identify the species/subspecies
  Usage: $0 run-dir
  --debug          Show debugging information
  --numcpus     1  Number of CPUs (has no effect on this script)
  --version
  "
}
