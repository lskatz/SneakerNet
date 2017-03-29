#!/usr/bin/env perl
# Use Kmers to guess if there is contamination

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;
use Bio::Kmer;

# http://perldoc.perl.org/perlop.html#Symbolic-Unary-Operators
# # +Inf and -Inf will be a binary complement of all zeros
use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help debug tempdir=s numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];
  
  my $report=kmerContaminationDetection($dir,$settings);

  # Write the output
  my $outfile="$dir/SneakerNet/forEmail/kmerHist.tsv";
  open(my $outFh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  print $outFh $report;
  close $outFh;
  logmsg "Kmer histograms were saved to $outfile";

  return 0;
}

sub kmerContaminationDetection{
  my($dir,$settings)=@_;
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my $i=0;

  my $report=join("\t",qw(sample hist firstPeak firstValley secondPeak secondValley...));
  logmsg $report;
  $report.="\n";
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(!ref($info));

    next if($$settings{debug} && rand() < 0.5);

    system("cat ".join(" ",@{$$info{fastq}})." > $$settings{tempdir}/$sample.fastq.gz");
    die "ERROR with cat into $$settings{tempdir}/$sample.fastq.gz" if $?;

    for my $fastq(@{$$info{fastq}}){
      logmsg "Counting kmers for $fastq";
      my $kmer=Bio::Kmer->new($fastq, {kmercounter=>"jellyfish",numcpus=>$$settings{numcpus}});
      my $hist=$kmer->histogram();

      # After kmer counting, find peaks and valleys
      my $peaksValleys=findThePeaksAndValleys($hist, 100, $settings);
      my $firstValley=$$peaksValleys{valleys}[0][0];
      my $lastValley =$$peaksValleys{valleys}[-1][0];
      if(!$firstValley || !$lastValley || $firstValley==$lastValley){
        logmsg "WARNING: no valleys detected in $fastq";
        $report.="$fastq\n";
        next;
      }

      # Look at the histogram between the first and last valley
      my @subHist=@$hist[$firstValley..$lastValley];
      my $histLine=$fastq;
      $histLine.="\t".sparkString(\@subHist,$settings);
      for(my $i=0;$i<@{$$peaksValleys{valleys}};$i++){
        $histLine.="\t".$$peaksValleys{peaks}[$i][0];
        $histLine.="\t".$$peaksValleys{valleys}[$i][0];
      }
      logmsg $histLine;
      $report.=$histLine."\n";

      if(++$i>3 && $$settings{debug}){
        last;
      }
    }
  }

  $report.="# The histogram ranges from the first valley to the last valley\n";
  $report.="# Some files have no discernable peaks and valleys and therefore might not have data listed\n";
  return $report;
}

sub findThePeaksAndValleys{
  my($hist, $delta, $settings)=@_;

  $settings||={};
  $$settings{gt}//=0;
  $$settings{maxzeros}//=3;

  my($min,$max)=(MAXINT,MININT);
  my($minPos,$maxPos)=(0,0);
  my @maxTab=();
  my @minTab=();

  my $lookForMax=1;

  my $numZeros=0; # If we see too many counts of zero, then exit.
  
  for(my $kmerCount=$$settings{gt}+1;$kmerCount<@$hist;$kmerCount++){
    my $countOfCounts=$$hist[$kmerCount];
    if($countOfCounts == 0){ 
      $numZeros++;
    }
    if($countOfCounts > $max){
      $max=$countOfCounts;
      $maxPos=$kmerCount;
    }
    if($countOfCounts < $min){
      $min=$countOfCounts;
      $minPos=$kmerCount;
    }

    if($lookForMax){
      if($countOfCounts < $max - $delta){
        push(@maxTab,[$maxPos,$max]);
        $min=$countOfCounts;
        $minPos=$kmerCount;
        $lookForMax=0;
      }
    }
    else{
      if($countOfCounts > $min + $delta){
        push(@minTab,[$minPos,$min]);
        $max=$countOfCounts;
        $maxPos=$kmerCount;
        $lookForMax=1;
      }
    }

    last if($numZeros > $$settings{maxzeros});
  }

  return {peaks=>\@maxTab, valleys=>\@minTab};
}

sub sparkString{
  my($arr,$settings)=@_;
  $settings||={};
  # Don't allow small numbers
  return "" if(ref($arr) ne "ARRAY" || @$arr < 2);
  # Make the string for spark
  my $histString=join(" ",@$arr);
  my $sparkString=`echo $histString | spark`;
  die "ERROR with spark" if $?;

  return $sparkString;
}


sub usage{
  "Guesses if there is contamination in a miseq run using kmer magic
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --debug       Just run three random samples
  "
}

