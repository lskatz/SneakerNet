#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

diag `baseBalance.pl --check-dependencies 2>&1`;
if($?){
  plan skip_all=>"Dependencies not met for baseBalance.pl";
}
plan tests=>2;

is system("baseBalance.pl --numcpus 1 --force $run >/dev/null 2>&1"), 0, "Running baseBalance.pl";

my $basebalance = "$run/SneakerNet/forEmail/basebalance.tsv";

# Double check results
subtest "Expected base balance results from $basebalance" => sub {
  plan tests => 20;
  my %AT = (
    "Philadelphia_CDC_1.fastq.gz" => 0.94,
    "Philadelphia_CDC_2.fastq.gz" => 0.94,
    "2010EL-1786_1.fastq.gz"      => 1.04,
    "2010EL-1786_2.fastq.gz"      => 1.01,
    "FA1090_1.fastq.gz"           => 0.91,
    "FA1090_2.fastq.gz"           => 0.91,
    "contaminated_1.fastq.gz"     => 1.15,
    "contaminated_2.fastq.gz"     => 1.15,
    "LT2_1.fastq.gz"              => 0.95,
    "LT2_2.fastq.gz"              => 0.93,
  );
  my %GC = (
    "Philadelphia_CDC_1.fastq.gz" => 1.12,
    "Philadelphia_CDC_2.fastq.gz" => 1.10,
    "2010EL-1786_1.fastq.gz"      => 0.99,
    "2010EL-1786_2.fastq.gz"      => 0.97,
    "FA1090_1.fastq.gz"           => 1.06,
    "FA1090_2.fastq.gz"           => 1.06,
    "contaminated_1.fastq.gz"     => 1.07,
    "contaminated_2.fastq.gz"     => 1.09,
    "LT2_1.fastq.gz"              => 0.97,
    "LT2_2.fastq.gz"              => 0.96,
  );

  open(my $fh, $basebalance) or die "ERROR reading base balance file $basebalance: $!";
  while(<$fh>){
    chomp;
    my ($file, $A, $T, $C, $G, $N, $AT, $GC)
        = split(/\t/, $_);
    
    next if(!$AT || !$AT{$file});

    is sprintf("%0.1f",$AT), sprintf("%0.1f", $AT{$file}), "AT percentage for $file";
    is sprintf("%0.1f",$GC), sprintf("%0.1f", $GC{$file}), "GC percentage for $file";
    #ok $AT==$AT{$file}, "AT percentage for $file";
    #ok $GC==$GC{$file}, "GC percentage for $file";
  }
  close $fh;
};
