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


subtest 'staramr' => sub{
  diag `sn_staramr.pl --check-dependencies 2>&1`;
  if($?){
    plan 'skip_all' => "sn_staramr.pl dependency check failed";
  }

  # Run staramr
  open(my $fh, "sn_staramr.pl $run --numcpus 2 --force 2>&1 | ") or BAIL_OUT("ERROR: could not run staramr on $run: $!");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  my $exit_code = $? << 8;
  is($exit_code, 0, "Running sn_staramr.pl");

  # Spot check the output
  subtest 'results' => sub{
    my $results = readTsv("$run/SneakerNet/forEmail/staramr.tsv");
    #is($$results{LT2}{"Predicted Phenotype"}, "Sensitive", "LT2 phenotype");
    isnt($$results{"2010EL-1786"}{"Predicted Phenotype"}, "Sensitive", "2010EL-1786 phenotype is not sensitive");
    ok($$results{"2010EL-1786"}{"Predicted Phenotype"} =~ /streptomycin|kanamycin|trimethoprim|chloramphenicol|sulfisoxazole"/i, "2010EL-1786 phenotype is at least streptomycin, kanamycin, trimethoprim, chloramphenicol, and/or sulfisoxazole");
  };
};
    
