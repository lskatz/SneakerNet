#!/usr/bin/env perl
# Use Kmers to guess if there is contamination

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;
use List::Util qw/min max/;

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
  GetOptions($settings,qw(help force debug tempdir=s numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/kmerHistogram";
  
  my $numFastq=kmerContaminationDetection($dir,$settings);
  
  # Write the output
  my $outfile="$dir/SneakerNet/forEmail/kmerHist.tsv";
  open(my $outFh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  print $outFh join("\t",qw(File numPeaks finalDelta hist firstPeak firstValley secondPeak secondValley...))."\n";
  # Collect all samples
  for my $fastqKmerCoverageFile (glob("$dir/SneakerNet/kmerHistogram/*/*.tsv")){
    open(my $fh, "<", $fastqKmerCoverageFile) or die "ERROR: could not read $fastqKmerCoverageFile: $!";
    my @line=<$fh>;
    close $fh;
    print $outFh @line;
  }
  print $outFh "# The histogram ranges from the first valley to the last valley\n";
  print $outFh "# Some files have no discernable peaks and valleys and therefore might not have data listed\n";
  close $outFh;

  logmsg "Kmer histograms were saved to $outfile";

  return 0;
}

sub kmerContaminationDetection{
  my($dir,$settings)=@_;
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  my $i=0;

  while(my($sample,$info)=each(%$sampleInfo)){
    next if(!ref($info));

    next if($$settings{debug} && rand() < 0.5);

    mkdir "$dir/SneakerNet/kmerHistogram/$sample";
    logmsg "$dir/SneakerNet/kmerHistogram/$sample";
    for my $fastq(@{$$info{fastq}}){
      logmsg "Counting kmers for $fastq";
      my $histGraph="$dir/SneakerNet/kmerHistogram/$sample/".basename($fastq).".graph.tsv";
      my $histTable="$dir/SneakerNet/kmerHistogram/$sample/".basename($fastq).".tsv";
      next if(!$$settings{force} && -e $histGraph);

      my $kmer=Bio::Kmer->new($fastq, {kmercounter=>"jellyfish",numcpus=>$$settings{numcpus}});
      my $hist=$kmer->histogram();

      # Write the histogram to disk
      open(my $histFh, ">", $histTable) or die "ERROR: could not write to $histTable: $!";
      for(my $i=0; $i<@$hist;$i++){
        print $histFh join("\t",$i,$$hist[$i])."\n";
      }
      close $histFh;

      # After kmer counting, find peaks and valleys
      my $peaksValleys=findThePeaksAndValleys($hist, 100, $settings);
      my $firstValley=$$peaksValleys{valleys}[0][0];
      my $lastValley =$$peaksValleys{valleys}[-1][0];
      if(!$firstValley || !$lastValley || $firstValley==$lastValley){
        logmsg "WARNING: no valleys detected in $fastq";
        next;
      }

      # Look at the histogram between the first and last valley
      my @subHist=@$hist[$firstValley..$lastValley];
      my $histLine=join("\t",basename($fastq),$$peaksValleys{numPeaks}, $$peaksValleys{delta});
      $histLine.="\t".sparkString(\@subHist,$settings);
      for(my $i=0;$i<@{$$peaksValleys{valleys}};$i++){
        $histLine.="\t".$$peaksValleys{peaks}[$i][0];
        $histLine.="\t".$$peaksValleys{valleys}[$i][0];
      }
      logmsg $histLine;
      open(my $fh, ">", $histGraph) or die "ERROR: could not write to $histGraph: $!";
      print $fh $histLine."\n";
      close $fh;

      # Do some cleanup
      $kmer->close;

      if(++$i>3 && $$settings{debug}){
        last;
      }
    }
  }

  return $i;

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
  # TODO exit the loop on a smart kmer coverage, like 2x expected genome coverage
  
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

  # If the delta is still a high enough number, but
  # there aren't at least two valleys, reiterate with
  # a smaller delta.
  if(!$minTab[0][0] || $minTab[0][0] == $minTab[-1][0]){
    if($delta > 10){
      logmsg "WARNING: I did not detect any peaks with delta=$delta.";
      $delta = int($delta / 2);
      logmsg "Reducing delta to $delta";
      return findThePeaksAndValleys($hist,$delta,$settings);
    } else {
      logmsg "WARNING: despite reducing the delta to $delta, no peaks were detected";
    }
  }

  return {peaks=>\@maxTab, valleys=>\@minTab, delta=>$delta, numPeaks=>scalar(@maxTab), numValleys=>scalar(@minTab)};
}

sub sparkString{
  my($values,$settings)=@_;

  my @ticks = qw/ ▁ ▂ ▃ ▄ ▅ ▆ ▇ █ /;

  my $min   = min @$values;
  my $max   = max @$values;
  my $range = $max - $min;
  my $size  = ($range << 8) / (scalar @ticks - 1);

  $size = 1 if $size < 1;

  my $sparkStr="";
  for my $v ( @$values ) {
    my $v_size = (($v - $min) << 8) / $size;

    $sparkStr.=$ticks[$v_size];
  }

  return $sparkStr;
}



sub usage{
  "Guesses if there is contamination in a miseq run using kmer magic
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --debug       Just run three random samples
  --force       Overwrite all results
  "
}

