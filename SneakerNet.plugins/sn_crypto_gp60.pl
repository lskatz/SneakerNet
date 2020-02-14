#!/usr/bin/env perl
# Runs GP60 on a genome assembly

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.0";
our $CITATION= "Gp60 plugin by Lee Katz. Uses gp60 by Alyssa Kelley.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      blastn    => 'blastn -version | head -n 1',
      rm        => 'rm --version | head -n 1',
      'countGP60repeats.pl' => 'echo GP60 counter unknown version',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/gp60";
  mkdir "$dir/SneakerNet/forEmail";
  my $reportArr=gp60ForAll($dir,$settings);
  my $reportFile = "$dir/SneakerNet/forEmail/gp60.tsv";

  # Make a report for email
  die;

  recordProperties($dir,{version=>$VERSION, table=>$reportFile});

  return 0;
}

sub gp60ForAll{
  my($dir, $settings) = @_;

  my @queueBuffer = ();

  # Find information about each genome
  logmsg "Reading sample tsv at $dir/samples.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    
    my $assembly=(glob("$dir/SneakerNet/assemblies/$sample/$sample.*.fasta"))[0];
    if(!$assembly){
      logmsg "ERROR: no assembly was found for $sample. Will not run analysis on this one sample.";
      next;
    }

    push(@queueBuffer, {assembly=>$assembly, sample=>$sample});
  }

  # Kick off the threads with the array buffer
  my $Q=Thread::Queue->new(@queueBuffer);
  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i]=threads->new(\&gp60Worker,$Q,$dir,$settings);
    # also send a terminating signal
    $Q->enqueue(undef);
  }

  my @reportArr=();
  for(@thr){
    my $reportList=$_->join;
    push(@reportArr, @$reportList);
  }

  return \@reportArr;
}

sub usage{
  print "Run gp60 analysis on a set of assemblies
  Usage: $0 MiSeq_run_dir
  --numcpus 1
";
  exit(0);
}

