#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;
use Scalar::Util /looks_like_number/;

use threads;

use Test::More tests => 2;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

my $kraken = `which kraken 2>/dev/null`; 
chomp($kraken);

subtest 'kraken results' => sub{
  if(!$kraken){
    plan skip_all => 'kraken executable not found';
  }

  plan tests => 37;

  diag `sn_detectContamination-kraken.pl --numcpus 2 --force $run 2>&1`;
  is $?, 0, "Running Kraken plugin";

  my $spreadsheet = "$run/SneakerNet/forEmail/kraken.tsv";

  diag "Reading from $spreadsheet";

  # Based on Kalamari v3.7
  my %expect = (
    contaminated => {
      NAME                               => "contaminated",
      LABELED_TAXON                      => "Legionella",
      BEST_GUESS                         => "Legionella pneumophila",
      PERCENTAGE_OF_GENOME_IS_BEST_GUESS => 75.24,
      MAJOR_CONTAMINANT                  => "Neisseria gonorrhoeae",
      PERCENTAGE_CONTAMINANT             => 16.97,
    },
    FA1090       => {
      NAME                               => "FA1090",
      LABELED_TAXON                      => "Neisseria",
      BEST_GUESS                         => "Neisseria gonorrhoeae",
      PERCENTAGE_OF_GENOME_IS_BEST_GUESS => 69.25,
      MAJOR_CONTAMINANT                  => ".",
      PERCENTAGE_CONTAMINANT             => 0,
    },
    Philadelphia_CDC => {
      NAME                               => "Philadelphia_CDC",
      LABELED_TAXON                      => "Legionella",
      BEST_GUESS                         => "Legionella pneumophila",
      PERCENTAGE_OF_GENOME_IS_BEST_GUESS => 99.77,
      MAJOR_CONTAMINANT                  => ".",
      PERCENTAGE_CONTAMINANT             => 0,
    },
    '2010EL-1786'    => {
      NAME                               => "contaminated",
      LABELED_TAXON                      => "Vibrio",
      BEST_GUESS                         => "Vibrio cholerae",
      PERCENTAGE_OF_GENOME_IS_BEST_GUESS => 82.47,
      MAJOR_CONTAMINANT                  => ".",
      PERCENTAGE_CONTAMINANT             => 0,
    },
  );

  open(KRAKEN, '<', $spreadsheet) or BAIL_OUT("ERROR: could not read from $spreadsheet: $!");
  my $header = <KRAKEN>;
  chomp($header);
  my @header = split(/\t/, $header);
  while(<KRAKEN>){
    next if(/^#/);

    chomp;
    my @F = split /\t/;
    my %F;
    @F{@header} = @F;

    my $name = $F{NAME};

    # Check some basic labels
    for my $h(qw(LABELED_TAXON BEST GUESS)){
      my $got = $F{$h};
      my $expected = $expect{$name}{$h};
      is($got, $expected, "$name $h");
    }

    # Check that some fields are numbers > 0
    for my $h(qw(PERCENTAGE_OF_GENOME_IS_BEST_GUESS PERCENTAGE_CONTAMINANT)){
      my $got = $F{$h};
      my $expected = $expect{$name}{$h};
      is(Scalar::Util::looks_like_number($got), 1, "Make sure $h is a number");
      cmp_ok($got, '>=', 0, "$name $h (check for number > 0)");
      cmp_ok($got, '<=', 100, "$name $h (check for number <= 100)");
    }
  }
  close KRAKEN;
};

