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
  pass("Skipping sars-cov-2 read metrics test in CI environment");
  exit 0;
}

my $numcpus = 2;
if($ENV{DEBUG}){
  $numcpus = 12;
  note "DEBUG: using $numcpus cpus";
}

subtest 'addReadMetrics' => sub {

  diag `addReadMetrics.pl --check-dependencies 2>&1`;
  if($?){
    plan 'skip_all' => "Plugin addReadMetrics.pl dependencies not met";
  } else {
    plan tests=>3;
  }

  is $?, 0, "addReadMetrics.pl dependency check";
  my $readMetricsLog = `addReadMetrics.pl --numcpus $numcpus --force $run 2>&1`;
  is $?, 0, "Adding read metrics";
  if($?){
    note $readMetricsLog;
  }

  subtest 'expected coverage plus or minus 10x' => sub{
    my $tolerance = 10;
    my %expected = (
      'SRR11826835' => 1193,
      'SRR12530737' => 1766,
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

