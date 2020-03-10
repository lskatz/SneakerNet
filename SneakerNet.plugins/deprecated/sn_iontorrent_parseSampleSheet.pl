#!/usr/bin/env perl
BEGIN{die "DEPRECATED: this script is no longer needed";}
# Reformats a sample sheet

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use Cwd qw/realpath/;
use File::Temp;
use File::Copy qw/cp mv/;
use File::Find;
use FindBin;
use Config::Simple;
use JSON ();
use Encode qw/encode decode/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig command logmsg/;

our $VERSION = "2.0";
our $CITATION= "Ion torrent sample sheet parsing by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(citation version help numcpus=i debug tempdir=s force)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }
  if($$settings{citation}){
    print $CITATION."\n";
    return 0;
  }

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  # We have to know if the argument is a directory or file, and so
  # we need to resolve the link.
  if(-l $ARGV[0]){
    $ARGV[0] = readlink($ARGV[0]);
  }
  my($dir) = @ARGV;
  my $outfile = "$dir/samples.tsv";
  if(-e $outfile){
    if(!$$settings{force}){
      die "$outfile already found. Exiting.";
    }
  }

  my $runInfo = iontorrentRunInfo($dir, $settings);

  writeSamplesTsv($runInfo, $outfile, $settings);

  recordProperties($ARGV[0],{version=>$VERSION, samples=>$outfile});

  return 0;
}

sub iontorrentRunInfo{
  my($dir, $settings) = @_;

  my $jsonIn = "$dir/planned_run.json";
  if(! -e $jsonIn){
    die "ERROR: cannot find $jsonIn";
  }

  my $json=JSON->new;
  $json->utf8;           # If we only expect characters 0..255. Makes it fast.
  $json->allow_nonref;   # can convert a non-reference into its corresponding string
  $json->allow_blessed;  # encode method will not barf when it encounters a blessed reference
  $json->pretty;         # enables indent, space_before and space_after

  my $jsonStr = "";
  open(my $jsonFh, "<", $jsonIn) or die "ERROR reading $jsonIn: $!";
  {
    local $/ = undef;
    $jsonStr = <$jsonFh>;
  }
  close $jsonFh;

  # Need to check for valid utf8 or not
  eval{ my $strCopy = $jsonStr; decode('utf8', $strCopy, Encode::FB_CROAK) }
    or die "ERROR: $dir/planned_run.json yielded non-utf8 characters\nContents shown below:\n$jsonStr\n";

  my $runInfo = $json->decode($jsonStr);
  return $runInfo;
}

sub writeSamplesTsv{
  my($runInfo, $outfile, $settings) = @_;

  die "TODO";
  open(my $fh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";

  my %barcode;
  my $barcodedSamples = $$runInfo{objects}[0]{barcodedSamples};
  my @sample = sort keys(%$barcodedSamples);
  for my $sample(@sample){
    while(my($barcode, $barcodeHash) = each(%{ $$barcodedSamples{$sample}{barcodeSampleInfo} })){
      print "$barcode\t$barcodeHash\n";
    }
  }
  die;

=cut
  open(my $outFh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  for my $barcode(sort keys(%$samples)){
    my $taxon = $$samples{$barcode}{organism};
    $taxon=~s/\S+.*$//; # remove everything after whitespace
    print $outFh join("\t",
      $$samples{$barcode}{'sample name'},
      "taxon=$taxon;route=calcengine",
      $$samples{$barcode}{fastq},
    ) . "\n";
  }
  close $outFh;

=cut
  return $outfile;
}

sub usage{
  print "Reformat a SampleSheet.csv file in a SneakerNet run
  Usage: $0 run-dir
  
  Because this is a SneakerNet plugin, the first positional
  argument must be the SneakerNet run. This script
  will not accept the spreadsheet as the first
  argument.

  OUTPUT: this script will print the new samplesheet to 
  stdout.
";
  exit(0);
}
