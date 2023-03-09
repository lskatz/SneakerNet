#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 3;
use Data::Dumper;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readConfig/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

# This is a difficult test because we don't really know
# if the email box received an email.
subtest 'sn_immediateStatus.pl' => sub{

  if($ENV{CI}){
    plan 'skip_all' => "Detected CI environment. Skipping immediate status test.";
  }

  diag `sn_immediateStatus.pl --check-dependencies 2>&1`;
  if($?){
    plan skip_all=>"Dependencies not met for sn_immediateStatus.pl";
  }

  system("sn_immediateStatus.pl $run");
  my $exit_code = $? >> 8;
  is($exit_code, 0, "Exit code");
};

subtest 'sn_report.pl' => sub{
  diag `sn_report.pl --check-dependencies 2>&1`;
  if($?){
    plan skip_all=>"Dependencies not met for sn_report.pl";
  }

  diag `sn_report.pl $run 2>&1`;
  my $exit_code = $? >> 8;
  is($exit_code, 0, "Exit code");
};

subtest 'emailWhoever.pl' => sub{
  diag `emailWhoever.pl --check-dependencies 2>&1`;
  if($?){
    plan skip_all=>"Dependencies not met for emailWhoever.pl";
  }
  system("emailWhoever.pl $run");
  my $exit_code = $? >> 8;
  is($exit_code, 0, "Exit code");
};

