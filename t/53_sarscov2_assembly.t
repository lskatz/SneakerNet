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

if($ENV{CI}){
  pass("Skipping sars-cov-2 assembly test in CI environment");
  exit 0;
}

subtest 'assembly' => sub {

  # Test for the older version of bgzip which we can't handle
  diag `which bgzip tabix`;
  diag `bgzip --help 2>&1`;
  if($?){
    plan 'skip_all' => "bgzip is an older version.";
  }

  diag `sn_sarscov2_assembleAll.pl --check-dependencies 2>&1`;
  if($?){
    plan 'skip_all' => "Plugin sn_sarscov2_assembleAll.pl dependencies not met";
  } else {
    plan tests=>3;
  }
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
        altCdsCount   => 9.00,
        refCdsCount   => 9.00,
        percentNs     => 0.01,
      },
      'SRR12530737.bowtie2.bcftools.fasta' => {
        altCdsCount   => 9.00,
        refCdsCount   => 9.00,
        percentNs     => 0.01,
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
