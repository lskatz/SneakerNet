#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readTsv/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-20-001-sarscov2";

subtest 'assembly' => sub {

  plan tests => 3;
  diag `which bgzip tabix`;
  diag `sn_sarscov2_assembleAll.pl --check-dependencies 2>&1`;
  is $?, 0, "sn_sarscov2_assembleAll.pl dependency check";

  open(my $fh, "sn_sarscov2_assembleAll.pl --force --numcpus 12 $run 2>&1 | ");
  while(<$fh>){
    diag $_;
  }
  close $fh;
  is($?, 0, "Ran the assembly");

  # File                                GC     N50    altCdsCount  assemblyScore  avgContigLength  expectedCdsPercentage  expectedGenomeLength  genomeLength  kmer21  longestContig  minContigLength  numContigs  percentNs  refCdsCount
  # SRR11826835.bowtie2.bcftools.fasta  0%     29903  10           10.306         29903            1.11                   0                     29903         0.0000  29903          500              1           1.00       9
  # SRR12530737.bowtie2.bcftools.fasta  14.7%  29903  9            10.306         29903            1.00                   0                     29903         0.0033  29903          500              1           0.61       9
  
  subtest 'assembly metrics' => sub{
    my @metricsToTest = qw(altCdsCount refCdsCount percentNs);
    my $obs = readTsv("$run/SneakerNet/forEmail/assemblyMetrics.tsv");
    my %expected = (
      'SRR11826835.bowtie2.bcftools.fasta' => {
        altCdsCount   => 10.00,
        refCdsCount   => 9.00,
        percentNs     => 1.00,
      },
      'SRR12530737.bowtie2.bcftools.fasta' => {
        altCdsCount   => 9.00,
        refCdsCount   => 9.00,
        percentNs     => 0.61,
      },
    );

    plan tests => scalar(@metricsToTest) * scalar(keys(%expected));
    
    while(my($asm, $metrics) = each(%$obs)){
      for my $m(@metricsToTest){
        my $obs = sprintf("%0.2f", $$metrics{$m});
        my $exp = sprintf("%0.2f", $expected{$asm}{$m});
        is($obs, $exp, "$asm / $m");
      }
    }
  };
};
