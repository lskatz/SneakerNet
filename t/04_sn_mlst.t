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

is system("sn_mlst.pl --numcpus 1 --force $run >/dev/null 2>&1"), 0, "7-gene MLST";

my $mlstFile = "$run/SneakerNet/forEmail/mlst.tsv";

# Double check MLST calls
subtest "Expected MLST calls from $mlstFile" => sub {
  plan tests => 6;
  my %mlst_scheme = (
    "2010EL-1786.skesa.fasta"      => "vcholerae",
    "Philadelphia_CDC.skesa.fasta" => "-",         # Legionella is not yet defined
    "FA1090.skesa.fasta"           => "neisseria",
  );
  my %adk = (
    "2010EL-1786.skesa.fasta"      => 7,
    "Philadelphia_CDC.skesa.fasta" => "",          # Legionella is not yet defined
    "FA1090.skesa.fasta"           => 39,
  );
  open(my $fh, $mlstFile) or die "ERROR reading MLST file $mlstFile: $!";
  while(<$fh>){
    chomp;
    my ($file, $mlst_scheme, $ST, @locus) 
        = split(/\t/, $_);
    
    next if(!$mlst_scheme{$file}); # avoid header

    ok $mlst_scheme eq $mlst_scheme{$file}, "MLST scheme for $file";
    my $adk = "";
    for my $locus(@locus){
      if($locus =~ /adk.*?(\d+)/){
        $adk = $1;
      }
    }
    ok $adk eq $adk{$file}, "adk allele call for $file ($adk vs $adk{$file})";
  }
  close $fh;
};
