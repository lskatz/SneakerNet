#!/usr/bin/env perl
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

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help numcpus=i debug tempdir=s force)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  # We have to know if the argument is a directory or file, and so
  # we need to resolve the link.
  if(-l $ARGV[0]){
    $ARGV[0] = readlink($ARGV[0]);
  }

  my $outfile     = "";
  my $samplesheet = "";
  if(-d $ARGV[0]){
    $outfile = "$ARGV[0]/samples.tsv";
    $samplesheet = findSampleSheet($ARGV[0], $settings);
  } else {
    $samplesheet = $ARGV[0];
  }

  my $fastqs = prepareIonTorrentRun($ARGV[0],$settings);

  my $sampleHash = parseIonTorrentSampleSheet($samplesheet, $fastqs, $settings);

  writeSamplesTsv($sampleHash, $outfile, $settings);

  return 0;
}

sub findSampleSheet{
  my($dir, $settings)=@_;

  if(-e "$dir/SampleSheet.tsv"){
    return "$dir/SampleSheet.tsv";
  }

  die "ERROR could not find sample sheet in $dir";
}

sub prepareIonTorrentRun{
  my($dir, $settings) = @_;

  my @zip = glob("$dir/*.zip");

  for my $zip(@zip){
    my $basename = basename($zip);
    command("cd '$dir' && unzip -o '$basename'");
  }

  # Find fastq files, gzip them, collect them
  my %fastq;
  find({no_chdir=>1, wanted=>sub{
    my $path = $File::Find::name;
    return if(!-f $path);
    my $filename = basename($path);

    my $barcode = -1;
    if($path=~/_0*(\d+)\.f/){
      $barcode = $1;
    }

    my $newpath = "";

    if($path =~ /\.fastq$|\.fq$/){
      command("gzip -v '$path'");
      $newpath = "$dir/$filename.gz";
      mv("$path.gz", $newpath);
    }
    elsif($path =~ /\.fastq.gz$|\.fq.gz$/){
      $newpath = "$dir/$filename";
      mv($path, $newpath);
    }
    $fastq{$barcode} = $newpath;
  }}, $dir);

  #for my $zip(@zip){
    #unlink($zip);
  #}

  return \%fastq;
}

sub parseIonTorrentSampleSheet{
  my($samplesheet, $fastqs, $settings) = @_;

  my %sample;

  open(my $fh, "<", $samplesheet) or die "ERROR: could not read $samplesheet: $!";
  while(<$fh>){
    chomp;
    next if(!/Well position/);
    my @header = split(/\t/, lc $_);
    while(<$fh>){
      chomp;
      my %F;
      my @value = split(/\t/, $_);
      @F{@header} = @value;

      my $sampleIndex = $F{'p2 barcode'};
      last if(!$sampleIndex);

      $sample{$sampleIndex} = \%F;
      $sample{$sampleIndex}{fastq} = $$fastqs{$sampleIndex};
    }
  }
  close $fh;

  return \%sample;
}

sub writeSamplesTsv{
  my($samples, $outfile, $settings) = @_;

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

  return $outfile;
}

sub usage{
  "Reformat a SampleSheet.csv file in a SneakerNet run
  Usage: $0 run-dir|SampleSheet.csv
  
  Because this is a SneakerNet plugin, the first positional
  argument must be the SneakerNet run.  However, this script
  will also accept the spreadsheet itself as the first
  argument.

  OUTPUT: this script will print the new samplesheet to 
  stdout but if the first argument is a directory, it
  will place the file into the directory.
  "
}
