#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use List::Util qw/sum/;

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
my $exit_code = $? >> 8;
is $exit_code, 0, "Adding read metrics";

note $readMetricsLog;

# Double check that everything is 10x.  ish.
# The exception is that the Vibrio sample is only the
# 3 Mb chromosome and so the coverage calculation
# will be off.
note "Reading from $run/readMetrics.tsv";

subtest "Expected coverage" => sub {
  plan tests => 4;
  my %expected = (
    "FA1090"                      => [5, 5],
    "2010EL-1786"                 => [3.8, 3.8],
    "Philadelphia_CDC"            => [5, 5],
    "LT2"                         => [5, 5],
  );
  open(my $fh, "$run/readMetrics.tsv") or die "ERROR reading $run/readMetrics.tsv: $!";
  my $header = <$fh>;
  chomp($header);
  my @header = split(/\t/, $header);
  while(<$fh>){
    chomp;
    my %F;
    @F{@header} = split /\t/;

    my $sample = $F{Sample};

    if(!$expected{$sample}){ 
      note "Skipping sample $sample: not found in list of results to check against";
      next;
    }

    subtest "Coverage for $sample" => sub{
      plan tests => 3;
      my $expCoverage = sum(@{ $expected{$sample} });

      my $is_numerical = isnt($F{coverage}, '.', "Coverage is numerical");
      if($is_numerical){
        cmp_ok($F{coverage}, '>', $expCoverage - 1, "Coverage numerical test");
        cmp_ok($F{coverage}, '<', $expCoverage + 1, "Coverage numerical test");
      } else {
        fail("Coverage");
        fail("Coverage");
      }
    };
  }
  close $fh;
};

