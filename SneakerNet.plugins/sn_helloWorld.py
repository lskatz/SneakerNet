#!/usr/bin/env python3
# Hello World example

import tempfile
import argparse
import sys
import os

VERSION='1.0'
CITATION='Hello World example by Lee Katz'

def main(args):
  print(args);

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="SneakerNet Options", allow_abbrev=True)
  parser.add_argument('-v','--version',  action='version', version=VERSION);
  parser.add_argument('-c','--citation', action='store_true', help="Show the citation and exit");
  parser.add_argument('-t','--tempdir',  nargs=1, type=str, help="Set a temporary directory path. Default: tempfile.TemporaryDirectory()")
  args = parser.parse_args();
  if args.citation:
    print(CITATION)
    sys.exit(0)
    
  with tempfile.TemporaryDirectory() as tempdirname:
    if args.tempdir is None:
      args.tempdir = tempdirname
    else:
      args.tempdir = args.tempdir[0];
      print(args.tempdir)
      if os.path.isdir(args.tempdir):
        print("ERROR: "+args.tempdir+" already exists!")
        sys.exit(1)
      else:
        os.mkdir(args.tempdir)

    main(args)


"""
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

our $VERSION = "1.0";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version help force tempdir=s debug numcpus=i)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

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
  "Run Hello World as a plugin for SneakerNet.
  Creates a table in the 'forEmail' folder with which
  options have been specified.

  Usage: $0 MiSeq_run_dir
  --numcpus  1  Number of CPUs
  --force       Whether to overwrite results
  --debug       Print debugging information
  --version     Print the version and exit
  --help        This usage menu
  --tempdir  '' Specify the temporary directory
  "
}
"""

