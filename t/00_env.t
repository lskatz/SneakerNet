#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname/;

use Test::More tests => 3;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

subtest 'Perl modules' => sub{
  my @module = qw(threads Config::Simple);
  plan tests => scalar(@module);

  for my $module(@module){
    #system("perl -M$module -e 1");
    #is($?, 0, $module);
    eval "use $module;";
    ok(!$@, $module);
  }
};


$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";
my $run = "$RealBin/M00123-18-001-test";

# Set up the directory structure for the rest of the tests.
mkdir "$run/SneakerNet";
mkdir "$run/SneakerNet/forEmail";

# Make a contaminant genome
subtest 'Make a contaminated genome' => sub{
  plan tests => 5;
  system("zcat $run/Philadelphia_CDC_1.fastq.gz | gzip -c1 >  $run/contaminated_1.fastq.gz");
  is($?, 0, "Zcat Philadelphia_CDC_1.fastq.gz to contaminated R1");

  system("zcat $run/Philadelphia_CDC_2.fastq.gz | gzip -c1 >  $run/contaminated_2.fastq.gz");
  is($?, 0, "Zcat Philadelphia_CDC_2.fastq.gz to contaminated R2");

  my $numContaminantReads = 88888 * 2;
  system("zcat $run/FA1090_1.fastq.gz | head -n $numContaminantReads | gzip -c1 >> $run/contaminated_1.fastq.gz");
  is($?, 0, "Zcat $numContaminantReads lines of FA1090_1.fastq.gz to contaminated R1");
  system("zcat $run/FA1090_2.fastq.gz | head -n $numContaminantReads | gzip -c1 >> $run/contaminated_2.fastq.gz");
  is($?, 0, "Zcat $numContaminantReads lines of FA1090_2.fastq.gz to contaminated R2");

  # Copy the last line of the sample sheet but change the sample name
  open(SAMPLESHEET, '<', "$run/SampleSheet.csv") or die "ERROR reading $run/SampleSheet.csv: $!";
  my @header;
  my %F;
  my $need_to_add_contaminant = 1;
  while(<SAMPLESHEET>){
    chomp;
    my @F = split /,/;
    if(/Sample_ID/){
      @header = @F;
      next;
    }

    if(!@header){
      next;
    }

    @F{@header} = @F;

    if($F{Sample_ID} =~ /contaminated/i){
      $need_to_add_contaminant = 0;
    }
  }
  close SAMPLESHEET;

  # Avoid using the contaminant genome twice
  if($need_to_add_contaminant){
    my $newLine = "";
    open(SAMPLESHEET, '>>', "$run/SampleSheet.csv") or die "ERROR writing to $run/SampleSheet.csv: $!";
    for(@header){
      if(/Sample_ID/){
        $newLine .= "contaminated,";
      } else {
        $newLine .= $F{$_}.",";
      }
    }
    $newLine =~ s/,$//;
    print SAMPLESHEET $newLine."\n";
    close SAMPLESHEET;
  }

  my $numSamples = `grep -A 100 Sample_ID $run/SampleSheet.csv | tail -n +2 | wc -l`; 
  chomp($numSamples);
  is($numSamples, 5, "Added contaminant sample");
};

