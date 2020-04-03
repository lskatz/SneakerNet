#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;
use File::Copy qw/cp/;

use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

use Bio::SeqIO;
use IO::Compress::Gzip;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";

my $run = "$RealBin/M00123-18-003-metagenomics";

# Set up the directory structure for the rest of the tests.
mkdir $run;
mkdir "$run/SneakerNet";
mkdir "$run/SneakerNet/forEmail";

# Create a mock community
# Also create a "second" mock community that only has the first 200k reads of mock1
my $metagenomicsR1 = "$run/mock_1.fastq.gz";
my $metagenomicsR2 = "$run/mock_2.fastq.gz";
my $metagenomics2R1= "$run/mock2_1.fastq.gz";
my $metagenomics2R2= "$run/mock2_2.fastq.gz";
open(my $R1fh, ">", $metagenomicsR1) or die "ERROR: could not write to metagenomics R1 $metagenomicsR1: $!";
open(my $R2fh, ">", $metagenomicsR2) or die "ERROR: could not write to metagenomics R2 $metagenomicsR2: $!";
# R1 files
for my $fastq(glob("$RealBin/M00123-18-001-test/*_1.fastq.gz")){
  diag "cp $fastq";
  open(my $fastqFh, "<", $fastq) or die "ERROR: could not read fastq $fastq: $!";
  while(my $chunk=<$fastqFh>){
    print $R1fh $chunk;
  }
  close $fastqFh;
}

# R2 files
for my $fastq(glob("$RealBin/M00123-18-001-test/*_2.fastq.gz")){
  diag "cp $fastq";
  open(my $fastqFh, "<", $fastq) or die "ERROR: could not read fastq $fastq: $!";
  while(my $chunk=<$fastqFh>){
    print $R2fh $chunk;
  }
  close $fastqFh;
}
close $R1fh;
close $R2fh;

# Make the second mock community of only 200k reads
my $readCounter = 0;
open($R1fh, "zcat $metagenomicsR1 | ") or die "ERROR: could not read $metagenomicsR1: $!";
open(my $R1fh2, " | gzip -c > $metagenomics2R1") or die "ERROR: could not write to metagenomics R1 $metagenomics2R1: $!";
while(my $line = <$R1fh>){
  last if(++$readCounter > 800000);
  print $R1fh2 $line;
}
close $R1fh;
close $R1fh2;

$readCounter = 0;
open($R2fh, "zcat $metagenomicsR2 | ") or die "ERROR: could not read $metagenomicsR2: $!";
open(my $R2fh2, " | gzip -c > $metagenomics2R2") or die "ERROR: could not write to metagenomics R2 $metagenomics2R2: $!";
while(my $line = <$R2fh>){
  last if(++$readCounter > 800000);
  print $R2fh2 $line;
}
close $R2fh;
close $R2fh2;

# Add a commensal E. coli to mock1
my $commensalEcoliFasta = "$run/ecoli_MGY.fasta";
if( ! -e $commensalEcoliFasta ){
  diag `wget 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=CP019629&rettype=fasta' -O $commensalEcoliFasta 2>&1`;
  if($?){
    BAIL_OUT("ERROR downloading CP019629 commensal E. coli to $commensalEcoliFasta");
  }
}

# quick simulation of commensal Ecoli
# Write to memory first bc faster and bc small metagenome
my $step = 100;
my $seqlength = 250;
my $fragmentSize = 500;
my $qualLine = 'I' x $seqlength;

my $seqin=Bio::SeqIO->new(-file=>"$commensalEcoliFasta");
my $R1buffer;
my $R2buffer;
my $idCounter = 0;
while(my $seq=$seqin->next_seq){
  my $sequence = $seq->seq;
  my $length=length($sequence);
  my $calculatedLength = $length-$fragmentSize;
  for(my $i=0;$i<$calculatedLength;$i+=$step){
    my $R1seq = substr($sequence, $i, $seqlength);
    # revcom
    my $R2seq = substr($sequence, $i+$fragmentSize, $seqlength);
    $R2seq =~ tr/ATCGatcg/TAGCtagc/;
    $R2seq = reverse($R2seq);

    $idCounter++;
    $R1buffer.="\@CP019629MGY_$idCounter/1\n$R1seq\n+\n$qualLine\n";
    $R2buffer.="\@CP019629MGY_$idCounter/2\n$R2seq\n+\n$qualLine\n";

    if($idCounter % 10000 == 0){
      diag "$idCounter entries ...";
    }
  }
}

open(my $R1fh_2, ">", "$run/tmp_1.fastq") or die "ERROR: could not write to metagenomics R1 $run/tmp_1.fastq: $!";
open(my $R2fh_2, ">", "$run/tmp_2.fastq") or die "ERROR: could not write to metagenomics R2 $run/tmp_2.fastq: $!";
print $R1fh_2 $R1buffer;
print $R2fh_2 $R2buffer;
close $R1fh_2;
close $R2fh_2;
diag `gzip -c $run/tmp_1.fastq 2>&1 1>> $metagenomicsR1`;
diag `gzip -c $run/tmp_2.fastq 2>&1 1>> $metagenomicsR2`;
unlink("$run/tmp_1.fastq","$run/tmp_2.fastq");

open(my $snokFh, ">", "$run/snok.txt") or die "ERROR writing to snok.txt: $!";
print $snokFh "workflow = metagenomics\n";
close $snokFh;

open(my $samplesFh, ">", "$run/samples.tsv") or die "ERROR writing to samples.tsv: $!";
print $samplesFh "mock1\tTaxon=\tmock_1.fastq.gz;mock_2.fastq.gz\n";
print $samplesFh "mock2\tTaxon=\tmock2_1.fastq.gz;mock2_2.fastq.gz\n";
close $samplesFh;

subtest 'integrity of new run' => sub{
  my $r1LineCounter = 0;
  my $r2LineCounter = 0;

  # mock1
  open(my $fh1, "zcat $metagenomicsR1 |") or BAIL_OUT("ERROR: could not read $metagenomicsR1: $!");
  while(<$fh1>){
    $r1LineCounter++;
  }
  close $fh1;
  open(my $fh2, "zcat $metagenomicsR2 |") or BAIL_OUT("ERROR: could not read $metagenomicsR2: $!");
  while(<$fh2>){
    $r2LineCounter++;
  }
  close $fh2;
  
  my $expectedReads = 1634788;
  is($r1LineCounter, $expectedReads, "Num metagenomics R1 reads");
  is($r2LineCounter, $expectedReads, "Num metagenomics R2 reads");

  # mock2
  $r1LineCounter = 0;
  $r2LineCounter = 0;
  open($fh1, "zcat $metagenomics2R1 |") or BAIL_OUT("ERROR: could not read $metagenomics2R1: $!");
  while(<$fh1>){
    $r1LineCounter++;
  }
  close $fh1;
  open($fh2, "zcat $metagenomics2R2 |") or BAIL_OUT("ERROR: could not read $metagenomics2R2: $!");
  while(<$fh2>){
    $r2LineCounter++;
  }
  close $fh2;
  
  $expectedReads = 800000;
  is($r1LineCounter, $expectedReads, "Num metagenomics R1 reads");
  is($r2LineCounter, $expectedReads, "Num metagenomics R2 reads");

  open(my $snokFh, '<', "$run/snok.txt") or BAIL_OUT("ERROR: could not read $run/snok.txt: $!");
  while(<$snokFh>){
    chomp;
    my($key, $value) = split(/\s*=\s*/, $_);
    if($key =~ /workflow/){
      is($value, 'metagenomics', "Workflow");
      last;
    }
  }
  close $snokFh;

};

