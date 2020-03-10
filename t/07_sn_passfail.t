#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use Test::More tests => 3;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

is system("sn_passfail.pl --numcpus 1 --force $run >/dev/null 2>&1"), 0, "Running sn_passfail.pl";

my $passfail= "$run/SneakerNet/forEmail/passfail.tsv";

diag "Reading from $passfail";
diag `cat $passfail`;

# Double check results
subtest "Expected passfail results from $passfail" => sub {
  plan tests => 8;

  open(my $fh, $passfail) or die "ERROR reading passfail file at $passfail: $!";
  while(<$fh>){
    next if(/^#/);
    next if(/^Sample/);

    chomp;
    my ($file, $coverage, $quality)
        = split(/\t/, $_);
    
    is $coverage, 1, "Check if $file failed coverage";
    is $quality , 0, "Check if $file passed quality";
  }
  close $fh;
};
