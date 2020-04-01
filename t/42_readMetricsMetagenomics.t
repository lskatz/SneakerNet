#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-003-metagenomics";

my $readMetricsLog = `addReadMetrics.pl --force $run 2>&1`;
is $?, 0, "Adding read metrics";

note $readMetricsLog;

