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
use List::MoreUtils qw/uniq/;
use FindBin;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig logmsg samplesheetInfo command/;
use Email::Stuffer;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

# Global
my @nt=qw(A T C G N);
my @sampleHeader=("File",@nt,"A/T","C/G");

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug test numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  
  my $dir=$ARGV[0];

  my $out=baseBalanceAll($dir,$settings);
  
  command("cp -v $out $dir/SneakerNet/forEmail/");
  
  return 0;
}

sub baseBalanceAll{
  my($dir,$settings)=@_;

  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my $outdir="$dir/SneakerNet/basebalance";

  mkdir $outdir; die if $?;
  my @outfile;
  my $outfile="$outdir/basebalance.tsv";
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    my $balanceFile=baseBalance($s,$outdir,$settings);
    push(@outfile,$balanceFile);
  }
  
  open(OUTFILE,">",$outfile) or die "ERROR: could not open $outfile for writing: $!";
  print OUTFILE join("\t",@sampleHeader)."\n";
  close OUTFILE;
  command("grep -hv '^#' @outfile >> $outfile");

  return $outfile;
}

sub baseBalance{
  my($sHash,$outdir,$settings)=@_;

  logmsg $$sHash{sample_id};
  my $outfile="$outdir/$$sHash{sample_id}.tsv";;
  open(my $outFh,">",$outfile) or die "ERROR: could not open $outfile for writing: $!";

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
    $out.=$fastq."\t";
    $out.=$nt{$_}."\t" for(@nt);
    $out.=sprintf("%0.2f",($nt{A}/$nt{T}))."\t".sprintf("%0.2f",($nt{C}/$nt{G}))."\n";
    print $outFh $out;
  }
  close $outFh;
   
  return $outfile;
}


sub usage{
  "Runs the base-balance metric
  Usage: $0 runDir
  --numcpus 1
  "
}

