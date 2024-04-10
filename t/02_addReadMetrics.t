#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use Test::More tests => 3;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

diag `addReadMetrics.pl --check-dependencies 2>&1`;
if($?){
  BAIL_OUT("Plugin addReadMetrics.pl does not have all dependencies met");
}

my $readMetricsLog = `addReadMetrics.pl --force $run 2>&1`;
is $?, 0, "Adding read metrics";

note $readMetricsLog;

# Double check that everything is 10x.  ish.
# The exception is that the Vibrio sample is only the
# 3 Mb chromosome and so the coverage calculation
# will be off.
subtest "Expected coverage" => sub {
  plan tests => 8;
  my %expected = (
    "FA1090_1.fastq.gz"           => 5,
    "FA1090_2.fastq.gz"           => 5,
    "2010EL-1786_1.fastq.gz"      => 3.8,
    "2010EL-1786_2.fastq.gz"      => 3.8,
    "Philadelphia_CDC_1.fastq.gz" => 5,
    "Philadelphia_CDC_2.fastq.gz" => 5,
    "LT2_1.fastq.gz"              => 5,
    "LT2_2.fastq.gz"              => 5,
  );
  open(my $fh, "$run/readMetrics.tsv") or die "ERROR reading $run/readMetrics.tsv: $!";
  #system("wc -l $run/readMetrics.tsv");
  while(<$fh>){
    chomp;
    my ($file, $avgReadLength, $totalBases, $minReadLength, $maxReadLength, $avgQuality, $numReads, $coverage) 
        = split(/\t/, $_);

    $file = basename($file);
    
    # Skip header
    if($file =~ /File/i){
      next;
    }
    if(!$expected{$file}){ 
      note "Skipping fastq file $file: not found in list of results to check against";
      next;
    }
    subtest "Coverage for $file" => sub{
      plan tests => 3;
      my $is_numerical = isnt($coverage, '.', "Coverage is numerical");
      if($is_numerical){
        cmp_ok($coverage, '>', $expected{$file} - 1, "Coverage numerical test");
        cmp_ok($coverage, '<', $expected{$file} + 1, "Coverage numerical test");
      } else {
        fail("Coverage");
        fail("Coverage");
      }
    }
  }
  close $fh;
};

