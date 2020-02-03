#!/usr/bin/env perl
# Hello World example

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.1";
our $CITATION = "Hello world perl SneakerNet plugin by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation help force tempdir=s debug numcpus=i)) or die $!;
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
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

  if(! -d "$dir/SneakerNet"){
    mkdir "$dir/SneakerNet";
  }
  if(! -d "$dir/SneakerNet/forEmail"){
    mkdir "$dir/SneakerNet/forEmail";
  }

  my $outfile = makeTable($dir, $settings);
  logmsg "Table can be found at $outfile";

  recordProperties($dir,{version=>$VERSION, table=>$outfile});

  return 0;
}

sub makeTable{
  my($dir, $settings)=@_;

  my $outfile = "$dir/SneakerNet/forEmail/helloWorld.tsv";
  open(my $fh, ">", $outfile) or die "ERROR writing to $outfile: $!";
  print $fh join("\t",qw(key value))."\n";
  for my $key(qw(version help force tempdir debug numcpus)){
    $$settings{$key} ||= 0;
    print $fh join("\t", $key, $$settings{$key})."\n";
  }
  close $fh;

  return $outfile;
}
    

sub usage{
  print "Run Hello World as a plugin for SneakerNet.
  Creates a table in the 'forEmail' folder with which
  options have been specified.

  Usage: $0 MiSeq_run_dir
  --numcpus  1  Number of CPUs
  --force       Whether to overwrite results
  --debug       Print debugging information
  --version     Print the version and exit
  --help        This usage menu
  --tempdir  '' Specify the temporary directory
  ";
  exit 0;
}

