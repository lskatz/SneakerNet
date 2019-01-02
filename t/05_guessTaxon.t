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

my $kraken = `which kraken 2>/dev/null`; chomp($kraken);
if(! $kraken){
  diag "Kraken is not installed and so this whole unit test will be skipped";
  pass("kraken1");
  pass("kraken2");
  exit 0;
}

is system("guessTaxon.pl --numcpus 1 --force $run >/dev/null 2>&1"), 0, "Running guessTaxon.pl";

my $krakenFile = "$run/SneakerNet/forEmail/kraken.tsv";

# Double check results
subtest "Expected Kraken results from $krakenFile" => sub {
  plan tests => 3;
  my %percentage = (
    "FA1090"           => 50.26,
    "Philadelphia_CDC" => 99.53,
    "2010EL-1786"      => 95.01,
  );
  open(my $fh, $krakenFile) or die "ERROR reading Kraken file $krakenFile: $!";
  while(<$fh>){
    chomp;
    my ($name, $labeled, $guess, $percentage)
        = split(/\t/, $_);
    
    next if(!$percentage{$name} || !$percentage); # avoid header

    $percentage=~s/[^\d\.]//g; # only retain number characters

    ok $percentage > $percentage{$name} -1 && $percentage < $percentage{$name} +1, "Correct percentage of best guess for $name";
  }
  close $fh;
};
