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

my $skesa = `which skesa 2>/dev/null`; chomp($skesa);
if(! $skesa){
  diag "Skesa is not installed and so this whole unit test will be skipped";
  pass("assembly1");
  pass("assembly2");
  exit 0;
}


is system("assembleAll.pl --numcpus 1 --force $run >/dev/null 2>&1 "), 0, "Assembling all";

# Double check assembly metrics.
# Let the checks be loose though because of different
# versions of assemblers.
subtest "Expected assembly stats" => sub {
  plan tests => 6;
  my %genomeLength = (
    "2010EL-1786.skesa"      => 2934482,
    "Philadelphia_CDC.skesa" => 3287023,
    "FA1090.skesa"           => 1802838,
  );
  my %CDS = (
    "2010EL-1786.skesa"      => 2812,
    "Philadelphia_CDC.skesa" => 3223,
    "FA1090.skesa"           => 2122,
  );
  open(my $fh, "$run/SneakerNet/forEmail/assemblyMetrics.tsv") or die "ERROR reading $run/SneakerNet/forEmail/assemblyMetrics.tsv: $!";
  while(<$fh>){
    chomp;
    my ($file,$genomeLength,$CDS,$N50,$longestContig,$numContigs,$avgContigLength,$assemblyScore,$minContigLength,$expectedGenomeLength,$kmer21,$GC)
        = split(/\t/, $_);
    
    next if(!$genomeLength{$file}); # avoid header

    # Tolerance of 10k assembly length diff
    ok $genomeLength > $genomeLength{$file} - 10000 && $genomeLength < $genomeLength{$file} + 10000, "Genome length for $file";
    # tolerance of 50 CDS
    ok $CDS > $CDS{$file} - 50 && $CDS < $CDS{$file} + 50, "CDS count for $file";
  }
  close $fh;
};
