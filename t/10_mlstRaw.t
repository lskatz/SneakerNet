#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use threads;

use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readTsv/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";


subtest 'ColorID' => sub{
  if($ENV{CI}){
    plan 'skip_all' => "Detected CI environment. Skipping ColorID.";
  }

  diag `sn_detectContamination-mlst.pl --check-dependencies 2>&1`;
  if($?){
    diag "Dependencies not found. Will not test.";
    plan 'skip_all' => "sn_detectContamination-mlst.pl dependency check failed";
  }

  my $is_opened = open(my $fh, "sn_detectContamination-mlst.pl $run --numcpus 2 --force 2>&1 | ");
  ok(defined($is_opened), "Run the command sn_detectContamination-mlst.pl ($!)");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  is($?, 0, "Running sn_detectContamination-mlst.pl");

  subtest 'results' => sub{
    my $results = readTsv("$run/SneakerNet/forEmail/mlst-contamination-detection.tsv");
    is($$results{FA1090}{Scheme}, "neisseria", "Scheme for FA1090");
    is($$results{LT2}{Scheme}, "senterica", "Scheme for LT2");
    is($$results{LT2}{NumLociFound}, 7, "Loci found for LT2");
  };
};
    
