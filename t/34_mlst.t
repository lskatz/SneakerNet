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
  plan tests => 2;
  my $mlst = eval{fullPathToExec('mlst');};
  isnt($mlst, '', "MLST path: $mlst");
  open(my $fh, "sn_mlst.pl $run --numcpus 2 --force 2>&1 | ") or BAIL_OUT("ERROR: could not run mlst on $run: $!");
  while(my $msg = <$fh>){
    chomp($msg);
    diag $msg;
  }
  close $fh;
  is($?, 0, "Running sn_mlst.pl");
};
    
subtest 'wgMLST' => sub{
  pass("TODO FIND PATH to chewBBACA");
  pass("TODO run wgMLST");
};

