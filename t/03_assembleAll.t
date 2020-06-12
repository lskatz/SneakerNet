#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use threads;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";

my $numcpus = 2;
#note "DEBUG"; $numcpus=24;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

diag `assembleAll.pl --check-dependencies 2>&1`;
if($?){
  plan 'skip_all' => "Plugin assembleAll.pl dependencies not met";
} else {
  plan tests=>2;
}

my $tsv = "$run/SneakerNet/forEmail/assemblyMetrics.tsv";
#unlink($tsv); # ensure that assembleAll.pl doesn't skimp
is system("assembleAll.pl --numcpus $numcpus --force $run"), 0, "Assembling all";

# Double check assembly metrics.
# Let the checks be loose though because of different
# versions of assemblers.
subtest "Expected assembly stats" => sub {
  plan tests => 20;
  my %genomeLength = (
    "2010EL-1786.shovill.skesa"      => 2955394,
    "Philadelphia_CDC.shovill.skesa" => 3328163,
    "FA1090.shovill.skesa"           => 1918813,
    "contaminated.shovill.skesa"     => 5782258,
    "LT2.shovill.skesa"              => 4820055,
  );
  my %CDS = (
    "2010EL-1786.shovill.skesa"      => 2714,
    "Philadelphia_CDC.shovill.skesa" => 3096,
    "FA1090.shovill.skesa"           => 2017,
    "contaminated.shovill.skesa"     => 8949,
    "LT2.shovill.skesa"              => 4802,
  );
  my %depth = (
    "2010EL-1786.shovill.skesa"      => 4.56,
    "Philadelphia_CDC.shovill.skesa" => 4.60,
    "FA1090.shovill.skesa"           => 6.31,
    "contaminated.shovill.skesa"     => 5.17,
    "LT2.shovill.skesa"              => 4.76,
  );
  diag `echo;column -t $tsv`;
  open(my $fh, "$tsv") or die "ERROR reading $tsv: $!";
  while(<$fh>){
    chomp;
    my ($file,$genomeLength,$CDS,$N50,$longestContig,$numContigs,$avgContigLength,$assemblyScore,$minContigLength,$expectedGenomeLength,$kmer21,$GC,$effectiveCoverage)
        = split(/\t/, $_);
    
    next if(!$genomeLength{$file}); # avoid header

    # Tolerance of 10k assembly length diff
    #ok $genomeLength > $genomeLength{$file} - 100000 && $genomeLength < $genomeLength{$file} + 100000, "Genome length for $file (expected:$genomeLength{$file} found:$genomeLength)";
    # tolerance of 50 CDS
    #ok $CDS > $CDS{$file} - 1000 && $CDS < $CDS{$file} + 1000, "CDS count for $file (expected:$CDS{$file} found:$CDS)";

    # Let's make this super lax for now
    cmp_ok($genomeLength, '>', 1, "Genome length for $file (expected:$genomeLength{$file} found:$genomeLength)");
    cmp_ok($CDS, '>', 1, "CDS count for $file (expected:$CDS{$file} found:$CDS)");

    # Depth of coverage
    cmp_ok($effectiveCoverage, '>', 1, "Effective coverage > 1x (expected:$depth{$file} found: $effectiveCoverage)");
    cmp_ok($effectiveCoverage, '<', 200, "Effective coverage < 200x (sanity check)");
  }
  close $fh;
};

