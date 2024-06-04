#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;
use Scalar::Util qw/looks_like_number/;

use threads;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-003-metagenomics";

if($ENV{CI}){
  plan 'skip_all' => "Detected CI environment. Skipping metagenomics.";
}
plan tests=>2;

my $numcpus = 2;
if($ENV{DEBUG}){
  $numcpus = 12;
  note "DEBUG: using $numcpus cpus";
}

subtest 'kraken' => sub {
  diag `sn_kraken.pl --check-dependencies 2>&1`;
  if($?){
    plan 'skip_all' => "sn_kraken.pl dependency check failed";
  }
  # run kraken and print log messages as it goes
  open(my $fh, "sn_kraken.pl --numcpus $numcpus --force $run 2>&1 | ") or BAIL_OUT("ERROR: running Kraken plugin: $!");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  is $?, 0, "Running Kraken plugin";
};


subtest 'kraken results' => sub{
  pass("TODO future metagenomics kraken results plugin");
};

