#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use File::Path qw/remove_tree/;
use File::Copy qw/cp/;

use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $rawReadsDir = "$RealBin/M00123-18-001-test";
my $dir         = "$RealBin/M00123-18-001-asm";

subtest 'Create assembly-only dataset' => sub{
  remove_tree($dir);
  mkdir($dir);
  
  # Create the sample sheet
  my $opened = open(my $fh, '>', "$dir/samples.tsv");
  is($opened, 1, "Opened $dir/samples.tsv");

  # Grab all samples from the assembly folder of the raw sequences test folder
  for my $asm(glob("$rawReadsDir/SneakerNet/assemblies/*/*.fasta")){
    # Copy over the assembly file
    my $filename = basename($asm);
    my $target = "$dir/$filename";
    my $copied = cp($asm, $target);
    is($copied, 1, "Copying $asm => $target");

    # Add onto the sample sheet
    print $fh join("\t", basename($filename, ".skesa.fasta"), ".", $filename)."\n";
  }
  close $fh;

  # Create snok.txt
  my $opened2 = open(my $snokFh, '>', "$dir/snok.txt");
  is($opened, 1, "Opened $dir/snok.txt");
  print $snokFh "workflow = assembly\n";
  close $snokFh;
};

