#!/usr/bin/env perl
# Runs the NCBI "base-balance" metric which checks for whether
# there is approximately a 1 to 1 balance of As to Ts in the
# read set and a 1 to 1 balance between C and G.
# An offset of 1.3 is usually an indicator of problems.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/mv/;
use FindBin;

use threads;
use Thread::Queue;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig logmsg samplesheetInfo_tsv command/;
use List::MoreUtils qw/uniq/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

our $VERSION = "1.1";
our $CITATION = "Base Balance by Lee Katz";

# Global
my @nt=qw(A T C G N);
my @sampleHeader=("File",@nt,"A/T","C/G");

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation tempdir=s help inbox=s debug test force numcpus=i)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }
  if($$settings{citation}){
    print $CITATION."\n";
    return 0;
  }

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  
  my $dir=$ARGV[0];
  mkdir "$dir/forEmail";

  my $out=baseBalanceAll($dir,$settings);
  
  command("cp -v $out $dir/SneakerNet/forEmail/ >&2");

  recordProperties($dir,{version=>$VERSION, table=>"$dir/SneakerNet/forEmail/".basename($out)});
  
  return 0;
}

sub baseBalanceAll{
  my($dir,$settings)=@_;

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my $outdir="$dir/SneakerNet/basebalance";

  mkdir $outdir; die if $?;
  my @outfile;
  my $outfile="$outdir/basebalance.tsv";

  my @thr;
  my $Q=Thread::Queue->new();
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&baseBalanceWorker,$outdir,$Q,$settings);
  }
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    $Q->enqueue($s);
  }
  for(@thr){
    $Q->enqueue(undef);
  }
  for(@thr){
    my $balanceFile=$_->join();
    push(@outfile,@$balanceFile);
  }
  
  open(OUTFILE,">",$outfile) or die "ERROR: could not open $outfile for writing: $!";
  print OUTFILE join("\t",@sampleHeader)."\n";
  close OUTFILE;
  command("grep -hv '^#' @outfile >> $outfile");

  return $outfile;
}

sub baseBalanceWorker{
  my($outdir,$Q,$settings)=@_;
  
  my @bbFile=();
  while(defined(my $sHash=$Q->dequeue)){
    my $bbFile=baseBalance($sHash,$outdir,$settings);
    push(@bbFile,$bbFile);
  }
  return \@bbFile;
}

sub baseBalance{
  my($sHash,$outdir,$settings)=@_;

  logmsg $$sHash{sample_id};
  my $outfile="$outdir/$$sHash{sample_id}.tsv";;
  return $outfile if(-e $outfile && (-s $outfile > 0));

  open(my $outFh,">","$outfile.tmp") or die "ERROR: could not open $outfile.tmp for writing: $!";

  for my $fastq(@{ $$sHash{fastq} }){
    my $i=0;
    my %nt;
    $nt{$_}=1 for(@nt); # give a pseudocount just in case anything ends up being zero.
    open(my $fastqFh,"<",$fastq) or die "ERROR: could not open $fastq for reading: $!";
    while(<$fastqFh>){
      my $mod = ++$i % 4;
      if($mod==2){
        $_=uc($_);
        for my $nt(split(//)){
          $nt{$nt}++;
        }
      }
    }
    close $fastqFh;

    # header
    my $out='#'.join("\t",@sampleHeader)."\n";
    # values
    $out.=basename($fastq)."\t";
    $out.=$nt{$_}."\t" for(@nt);
    $out.=sprintf("%0.2f",($nt{A}/$nt{T}))."\t".sprintf("%0.2f",($nt{C}/$nt{G}))."\n";
    print $outFh $out;
  }
  close $outFh;

  mv("$outfile.tmp",$outfile) or die "ERROR: could not move $outfile.tmp to $outfile";
   
  return $outfile;
}


sub usage{
  print "Runs the base-balance metric
  Usage: $0 runDir
  --numcpus 1
  --version
  ";
  exit(0);
}

