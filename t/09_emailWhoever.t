#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 2;
use Data::Dumper;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readConfig/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

# This is a difficult test because we don't really know
# if the email box received an email.
subtest 'sn_immediateStatus.pl' => sub{
  system("sn_immediateStatus.pl $run");
  my $exit_code = $? >> 8;
  is($exit_code, 0, "Exit code");
};

subtest 'emailWhoever.pl' => sub{
  system("sn_immediateStatus.pl $run");
  my $exit_code = $? >> 8;
  is($exit_code, 0, "Exit code");
};

