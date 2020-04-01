#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use threads;

use Test::More tests => 2;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/fullPathToExec/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";


subtest 'Classic MLST' => sub{
  my $mlst = eval{fullPathToExec('mlst');};
  ok(defined($mlst), "mlst in path");

  diag `sn_mlst.pl --check-dependencies 2>&1`;
  if($?){
    fail("sn_mlst.pl dependencies check");
    return;
  } 

  open(my $fh, "sn_mlst.pl $run --numcpus 2 --force 2>&1 | ") or BAIL_OUT("ERROR: could not run mlst on $run: $!");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  is($?, 0, "Running sn_mlst.pl");
};
    
subtest 'wgMLST' => sub{
  diag `sn_mlst-wg.pl --check-dependencies 2>&1`;
  if($?){
    fail("sn_mlst-wg.pl dependency check failed");
    return;
  }

  open(my $fh, "sn_mlst-wg.pl $run --numcpus 2 --force 2>&1 | ") or BAIL_OUT("ERROR: could not run wgMLST on $run: $!");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  is($?, 0, "Running sn_mlst-wg.pl");
};

