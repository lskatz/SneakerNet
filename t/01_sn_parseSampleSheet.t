#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use Test::More tests => 2;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

# Set up the directory structure for the rest of the tests.
mkdir "$run/SneakerNet";
mkdir "$run/SneakerNet/forEmail";

system("sn_parseSampleSheet.pl --force $run >/dev/null 2>&1");
is $?, 0, "Parsing the sample sheet with sn_parseSampleSheet.pl";

