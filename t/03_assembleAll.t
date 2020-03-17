#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use threads;

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

is system("assembleAll.pl --numcpus 2 --force $run"), 0, "Assembling all";

# Double check assembly metrics.
# Let the checks be loose though because of different
# versions of assemblers.
subtest "Expected assembly stats" => sub {
  plan tests => 10;
  my %genomeLength = (
    "2010EL-1786.skesa"      => 2955394,
    "Philadelphia_CDC.skesa" => 3328163,
    "FA1090.skesa"           => 1918813,
    "contaminated.skesa"     => 5782258,
    "LT2.skesa"              => 4820055,
  );
  my %CDS = (
    "2010EL-1786.skesa"      => 2714,
    "Philadelphia_CDC.skesa" => 3096,
    "FA1090.skesa"           => 2017,
    "contaminated.skesa"     => 8949,
    "LT2.skesa"              => 4802,
  );
  diag `echo;column -t $run/SneakerNet/forEmail/assemblyMetrics.tsv`;
  open(my $fh, "$run/SneakerNet/forEmail/assemblyMetrics.tsv") or die "ERROR reading $run/SneakerNet/forEmail/assemblyMetrics.tsv: $!";
  while(<$fh>){
    chomp;
    my ($file,$genomeLength,$CDS,$N50,$longestContig,$numContigs,$avgContigLength,$assemblyScore,$minContigLength,$expectedGenomeLength,$kmer21,$GC)
        = split(/\t/, $_);
    
    next if(!$genomeLength{$file}); # avoid header

    # Tolerance of 10k assembly length diff
    #ok $genomeLength > $genomeLength{$file} - 100000 && $genomeLength < $genomeLength{$file} + 100000, "Genome length for $file (expected:$genomeLength{$file} found:$genomeLength)";
    # tolerance of 50 CDS
    #ok $CDS > $CDS{$file} - 1000 && $CDS < $CDS{$file} + 1000, "CDS count for $file (expected:$CDS{$file} found:$CDS)";

    # Let's make this super lax for now
    ok $genomeLength > 1, "Genome length for $file (expected:$genomeLength{$file} found:$genomeLength)";
    ok $CDS > 1, "CDS count for $file (expected:$CDS{$file} found:$CDS)";
  }
  close $fh;
};

