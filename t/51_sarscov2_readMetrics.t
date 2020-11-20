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

subtest 'addReadMetrics' => sub {

  plan tests => 3;
  diag `addReadMetrics.pl --check-dependencies 2>&1`;
  is $?, 0, "addReadMetrics.pl dependency check";
  my $readMetricsLog = `addReadMetrics.pl --numcpus 1 --force $run 2>&1`;
  is $?, 0, "Adding read metrics";
  if($?){
    note $readMetricsLog;
  }

  subtest 'expected coverage plus or minus 10x' => sub{
    my $tolerance = 10;
    my %expected = (
      'SRR11826835_1.fastq.gz' => 596,
      'SRR11826835_2.fastq.gz' => 596,
      'SRR12530737_1.fastq.gz' => 883,
      'SRR12530737_2.fastq.gz' => 883,
    );
    plan tests => scalar(keys(%expected)) * 2;
    my $obs = readTsv("$run/readMetrics.tsv");
    while(my($fastq, $metrics) = each(%$obs)){
      my $coverage = int($$metrics{coverage});
      my $expected = $expected{$fastq};
      cmp_ok($coverage, '<', $expected + $tolerance, "$fastq: $coverage < $expected+$tolerance");
      cmp_ok($coverage, '>', $expected - $tolerance, "$fastq: $coverage > $expected-$tolerance");
    }
  };

};

