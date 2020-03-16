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

is system("sn_passfail.pl --numcpus 1 --force $run >/dev/null 2>&1"), 0, "Running sn_passfail.pl";

my $passfail= "$run/SneakerNet/forEmail/passfail.tsv";

diag "Reading from $passfail";
diag `grep -v '^#' $passfail | column -t`;

# Double check results
subtest "Expected passfail results from $passfail" => sub {
  open(my $fh, $passfail) or die "ERROR reading passfail file at $passfail: $!";
  my $header = <$fh>;
  chomp($header);
  my @header = split(/\t/, $header);
  while(<$fh>){
    chomp;
    next if(/^#/);
    next if(/^Sample/);

    my @F = split(/\t/, $_);
    my %F;
    @F{@header} = @F;
    my $sample = $F{Sample};

    my $expectedContamination = 0;
    if($sample =~ /contamin/i){
      $expectedContamination = 1;
    }
    # Kraken is undefined if we never ran the test so...
    if(!-e "$run/SneakerNet/forEmail/kraken.tsv"){
      $expectedContamination = -1;
    }

    is $F{coverage}, 1, "Check if $sample failed coverage";
    is $F{quality} , 0, "Check if $sample passed quality";
    is $F{kraken}  , $expectedContamination, "Check if $sample passed read classification (kraken)";
  }
  close $fh;
};
