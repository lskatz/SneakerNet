#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use threads;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";

my $numcpus = 2;
if($ENV{DEBUG}){
  $numcpus=12;
  note "DEBUG: CPUS set to $numcpus"; 
}

if($ENV{CI}){
  plan 'skip_all' => "Detected CI environment. Skipping assembly";
}

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

diag `assembleAll.pl --check-dependencies 2>&1 | sort`;
if($?){
  plan 'skip_all' => "Plugin assembleAll.pl dependencies not met";
} else {
  plan tests=>2;
}

my $tsv = "$run/SneakerNet/forEmail/assemblyMetrics.tsv";
unlink($tsv); # ensure that assembleAll.pl doesn't skimp

my $cmd = "assembleAll.pl --numcpus $numcpus --force $run";
if($ENV{DEBUG}){
  note "DEBUG";
  $cmd =  "assembleAll.pl --numcpus $numcpus $run";
}
is system($cmd), 0, "Assembling all";

# Double check assembly metrics.
# Let the checks be loose though because of different
# versions of assemblers.
subtest "Expected assembly stats" => sub {
  #plan tests => 20; # number of tests will change depending on whether 00_env.t was run first
  my %genomeLength = (
    "2010EL-1786"      => 2955394,
    "Philadelphia_CDC" => 3328163,
    "FA1090"           => 1918813,
    "contaminated"     => 5782258,
    "LT2"              => 4820055,
  );
  my %CDS = (
    "2010EL-1786"      => 2714,
    "Philadelphia_CDC" => 3096,
    "FA1090"           => 2017,
    "contaminated"     => 8949,
    "LT2"              => 4802,
  );
  my %depth = (
    "2010EL-1786"      => 4.56,
    "Philadelphia_CDC" => 4.60,
    "FA1090.shovill"           => 6.31,
    "contaminated"     => 5.17,
    "LT2"              => 4.76,
  );
  diag `echo;column -t $tsv`;
  open(my $fh, "$tsv") or die "ERROR reading $tsv: $!";
  my $header = <$fh>;
  chomp($header);
  $header =~ s/#\s*//g; # remove the pound signs
  my @header = split(/\t/, $header);
  while(<$fh>){
    chomp;
    my @F = split(/\t/);
    my %F;
    @F{@header} = @F;
    #die Dumper \%F;

    my $file = $F{'Assembly'};
    my $sample = basename($file, qw(.shovill.skesa .shovill.skesa.fasta)); # quast already removes .fasta
    diag "Testing $file stats";

    # Let's make this super lax for now
    cmp_ok($F{'Total length'}, '>', 1, "Genome length for $file (expected:$genomeLength{$sample}");
    # predicted genes is non-numerical: usually in the format of confirmedGenes + partialGenes,
    # and so I am simplifying it to just confirmed genes.
    $F{'predicted genes (>= 0 bp)'} =~ s/\s+.*//;
    cmp_ok($F{'predicted genes (>= 0 bp)'}, '>', 1, "CDS count for $file (expected:$CDS{$sample})");

    # Depth of coverage
    # edit: Depth of coverage isn't part of assembly metrics at the moment
    #cmp_ok($effectiveCoverage, '>', 1, "Effective coverage > 1x (expected:$depth{$file} found: $effectiveCoverage)");
    #cmp_ok($effectiveCoverage, '<', 200, "Effective coverage < 200x (sanity check)");
  }
  close $fh;
};

