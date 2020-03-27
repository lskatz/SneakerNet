#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use File::Path qw/remove_tree/;
use File::Copy qw/cp/;

use Test::More tests => 2;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $dir         = "$RealBin/M00123-18-001-asm";

# Run the plugin
END{unlink("sn_assemblyWorkflow_init.pl.log");} # TODO just use tempdir
system("sn_assemblyWorkflow_init.pl --force $dir > sn_assemblyWorkflow_init.pl.log 2>&1");
my $exit_code = $? >> 8;
is($exit_code, 0, "Ran sn_assemblyWorkflow_init.pl");
if($exit_code){
  BAIL_OUT("ERROR running sn_assemblyWorkflow_init.pl and the log follows:\n".`cat sn_assemblyWorkflow_init.pl.log`);
}

# Check the output table, but the filenames are slightly different than the source
my $expected = `cut -f 2- M00123-18-001-test/SneakerNet/forEmail/assemblyMetrics.tsv | sort | md5sum`;
my $observed = `cut -f 2- M00123-18-001-asm//SneakerNet/forEmail/assemblyMetrics.tsv | sort | md5sum`;
is($observed, $expected, "Same assembly and gene prediction metrics as the original test data");

