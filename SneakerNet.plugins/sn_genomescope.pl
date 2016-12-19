#!/usr/bin/env perl
# Use Jellyfish/Genomescope to determine genome size

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use File::Copy qw/cp/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;

my $JELLYFISHDIR="/usr/local/bin";
my $GENOMESCOPEDIR="/opt/genomescope";
$ENV{PATH}="$ENV{PATH}:$JELLYFISHDIR:$GENOMESCOPEDIR";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help debug tempdir=s numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);
  $$settings{kmerlength}||=21;
  $$settings{readlength}||=150;

  my $dir=$ARGV[0];
  
  mkdir "$dir/SneakerNet/genomescope";
  predictGenomeSizesOnDir($dir,$settings);

  # make the report emailable 
  command("zip -v $dir/SneakerNet/forEmail/genomescope.zip $dir/SneakerNet/genomescope/* 1>&2");

  # Make spreadsheet
  my @field=("sample", "Heterozygosity", "Genome Haploid Length", "Genome Repeat Length", "Genome Unique Length", "Model Fit", "Read Error Rate");
  my $tsv=join("\t",@field)."\n";
  for my $genomescopeDir(glob("$dir/SneakerNet/genomescope/*")){
    my $summary="$genomescopeDir/summary.txt";
    open(SUMMARY,"<",$summary) or die "ERROR: could not open $summary: $!";
    #my $k=<SUMMARY>;
    #my $blank=<SUMMARY>;
    #my $header=<SUMMARY>;
    my %field;
    while(my $f=<SUMMARY>){
      chomp($f);
      my($property,$min,$max)=(split(/  +/,$f),"",0,0);
      $_=~s/[^\d\.%]//g for($min,$max); # remove nondigit, nondot characters
      $field{$property}="[$min - $max]";
    }
    $tsv.=basename(dirname($summary));
    for(my $i=1;$i<@field;$i++){
      my $header=$field[$i];
      # replace the percentage sign for some fields
      #if($header=~/Heterozygosity|Model Fit|Read Error Rate/){
      #  $field{$header}.='%';
      #}
      $tsv.="\t$field{$header}";
    }
    $tsv.="\n";
    close SUMMARY;
  }
  open(TSV,">","$dir/SneakerNet/genomescope.tsv") or die "ERROR: could not write to $dir/SneakerNet/genomescope.tsv: $!";
  print TSV $tsv;
  print TSV "NOTE: the genome haploid length is just an estimate of the genome size from raw reads.\n";
  print TSV "      The level of confidence is indicated by Model Fit.\n";
  print TSV "NOTE: the read error rate is an estimate of the reads' error rate using by kmer singletons.\n";
  close TSV;
  print $tsv;

  return 0;
}

sub predictGenomeSizesOnDir{
  my($dir,$settings)=@_;

  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    # Skip any samples without reads, ie, samples that are misnamed or not sequenced.
    # There is no way to predict how a sample is misnamed and so it does not fall under
    # this script's purview.
    if(!defined($$s{fastq}) || !@{ $$s{fastq} }){
      logmsg "WARNING: I could not find the reads for $sampleName. Skipping.";
      next;
    }
    
    my $histogram=runJellyfish($s,$settings);
    my $genomescopeDir=runGenomeScope($s,$histogram,$settings);

    # make results
    command("mv $genomescopeDir $dir/SneakerNet/genomescope/$sampleName");
    #logmsg "DEBUG";last;
  }
}

sub runJellyfish{
  my($sample,$settings)=@_;
  my $jfDir=tempdir("jellyfish.XXXX",DIR=>$$settings{tempdir});
  my $jfOut="$$settings{tempdir}/$$sample{sample_id}.jf";
  my $histoOut="$$settings{tempdir}/$$sample{sample_id}.histo";
  my $fastqFiles=join(" ",@{ $$sample{fastq} });

  # Check on the number of lines in the fastq file.
  # If it is less than 4, then there are zero reads.
  my @content=`zcat $fastqFiles | head -n 8888`; # pipe to head, to make it go fast
  my $numLines=scalar(@content);
  if($numLines < 888){
    logmsg "WARNING: $$sample{sample_id} has very few reads; I will fake a histogram for genome size";
    # make a fake histogram that Genomescope won't die on.
    open(HISTOUT,">", $histoOut) or die "ERROR: could not open $histoOut for writing: $!";
    for(1..50){
      my $count=50 - abs($_ - 25);
      $count=$count-25;
      $count=1 if($count < 1);
      print HISTOUT "$_ $count\n";
    }
    for(51..500){
      my $count=$_ % 3 + 1;
      print HISTOUT "$_ $count\n";
    }
    for(501..1000){
      my $count=$_ % 2 + 1;
      print HISTOUT "$_ $count\n";
    }
    close HISTOUT;
    return $histoOut;
  }

  # Make the histogram for genomescope
  logmsg "Kmer counting $$sample{sample_id} into $jfDir";
  command("zcat $fastqFiles | jellyfish count -C -m $$settings{kmerlength} -s 1000000000 -t $$settings{numcpus} -o $jfDir/JF /dev/stdin 1>&2");
  # Jellyfish outputs individual files with a given prefix.
  # However, they need to be merged into one file.
  my @jfFiles=glob("$jfDir/JF*");
  if(@jfFiles > 1){
    my $files=join(" ",@jfFiles);
    command("jellyfish merge -o $jfOut $files 1>&2");
  } else {
    cp($jfFiles[0],$jfOut);
  }
  # Run the histogram on the one output (merged) file.
  command("jellyfish histo -t $$settings{numcpus} $jfOut > $histoOut");
  return $histoOut;
}

sub runGenomeScope{
  my($sample,$histogram,$settings)=@_;
  
  my $genomeScopeOutdir=tempdir("genomescope.XXXX",DIR=>$$settings{tempdir});
  command("genomescope.R $histogram $$settings{kmerlength} $$settings{readlength} $genomeScopeOutdir 1>&2");
  return $genomeScopeOutdir;
}

sub usage{
  "Estimate the genome size using Jellyfish/Genomescope
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  "
}

