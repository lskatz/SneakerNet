#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-003-metagenomics";

diag `addReadMetrics.pl --check-dependencies 2>&1`;
if($?){
  BAIL_OUT("Plugin addReadMetrics.pl does not have all dependencies met");
}
plan tests=>1;

my $readMetricsLog = `addReadMetrics.pl --force $run 2>&1`;
is $?, 0, "Adding read metrics";

note $readMetricsLog;

