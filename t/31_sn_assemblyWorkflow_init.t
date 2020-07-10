#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use File::Path qw/remove_tree/;
use File::Copy qw/cp/;

use Test::More;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readTsv/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $dir         = "$RealBin/M00123-18-001-asm";

diag `sn_assemblyWorkflow_init.pl --check-dependencies 2>&1`;
if($?){
  plan 'skip_all' => "Dependencies not met for sn_assemblyWorkflow_init.pl";
}

plan tests => 2;

# Run the plugin
END{unlink("sn_assemblyWorkflow_init.pl.log");} # TODO just use tempdir
system("sn_assemblyWorkflow_init.pl --force $dir > sn_assemblyWorkflow_init.pl.log 2>&1");
my $exit_code = $? >> 8;
is($exit_code, 0, "Ran sn_assemblyWorkflow_init.pl");
if($exit_code){
  BAIL_OUT("ERROR running sn_assemblyWorkflow_init.pl and the log follows:\n".`cat sn_assemblyWorkflow_init.pl.log`);
}

if(!-e "$dir/SneakerNet/forEmail/assemblyMetrics.tsv"){
  plan 'skip_all' => "Assemblies not found. Skipping.";
}

# Check the output table, but the filenames are slightly different than the source
my $expected = readTsv("M00123-18-001-test/SneakerNet/forEmail/assemblyMetrics.tsv");
my $observed = readTsv("M00123-18-001-asm//SneakerNet/forEmail/assemblyMetrics.tsv");
$expected = sanitizeMetricsTable($expected);
$observed = sanitizeMetricsTable($observed);

is_deeply($observed, $expected, "Same assembly and gene prediction metrics as the original test data");

sub sanitizeMetricsTable{
  my($hash) = @_;

  my %h;

  while(my($filename, $metrics) = each(%$hash)){
    # Remove some fields that we don't exactly want to compare on
    delete($$metrics{effectiveCoverage});
    delete($$metrics{File});

    # Ignore extensions for the purposes of this test
    $filename =~ s/\..*//;

    # Copy over values instead of by reference
    my %newMetrics = %$metrics;
    $h{$filename} = \%newMetrics;
  }

  return \%h;

}

