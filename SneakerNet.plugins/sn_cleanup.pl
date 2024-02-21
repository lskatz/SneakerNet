#!/usr/bin/env perl
# Cleans up a SneakerNet folder

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use File::Path qw/rmtree/;
use List::Util qw/uniq/;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.1";
our $CITATION = "sn_cleanup.pl, by Lee Katz";

# For any warnings in the SN report
my $warningsMsg = "";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
    }, $settings,
  );

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{debug} ||= 0; # need this value to be explicitly defined

  my $dir=$ARGV[0];

  if(! -d "$dir/SneakerNet"){
    mkdir "$dir/SneakerNet";
  }
  if(! -d "$dir/SneakerNet/forEmail"){
    mkdir "$dir/SneakerNet/forEmail";
  }

  my %removalStats = cleanup($dir, $settings);
  my $numRemoved += $removalStats{num};
  my $sizeRemoved+= $removalStats{size};

  my $gigabytes = sprintf("%0.2f", $sizeRemoved / (1024 ** 3) );
  logmsg "Removed $numRemoved files to clear ${gigabytes}G";
  logmsg "NOTE: nothing was removed because --debug" if($$settings{debug});

  recordProperties($dir,{
    version=>$VERSION,
    warnings => $warningsMsg,
    debug => $$settings{debug}, 
    numRemoved => $numRemoved,
    sizeRemoved => "${gigabytes}G",
  });

  return 0;
}

sub cleanup{
  my($dir, $settings) = @_;

  my $numRemoved = 0;
  my $sizeRemoved= 0;

  # Which files do we remove?
  my @rmFiles = uniq glob(
    "
      $dir/SneakerNet/assemblies/*/shovill
      $dir/SneakerNet/assemblies/*/prodigal
    "
  );

  for my $file(@rmFiles){
    my %removalStats = cleanThis($file, $$settings{debug}, $settings);
    $numRemoved += $removalStats{num};
    $sizeRemoved+= $removalStats{size};
  }

  return (size=>$sizeRemoved, num=>$numRemoved);
}

sub cleanThis{
  my($file, $debug, $settings) = @_;

  my $numRemoved = 0;
  my $sizeRemoved = 0;

  if(-d $file){
    for my $f(glob("$file/*")){
      # But keep files ending in .log
      next if($f =~ /\.log$/);

      my %removalStats = cleanThis($f, $debug, $settings);
      $numRemoved += $removalStats{num};
      $sizeRemoved+= $removalStats{size};
    }

    # Now all intended files have been removed.
    # Remove the directory if it is empty.
    if(is_folder_empty($file)){
      logmsg "RMDIR: $file";
      $numRemoved++;
      my $size = -s $file;
      $sizeRemoved += $size;
      if(!$$settings{debug}){
        rmdir($file) or die "ERROR: could not rmdir $file: $!";
      }
    } else {
      logmsg "NOTE: $file/ not empty; not removing.";
    }
  }
  else {
    logmsg "UNLINK: $file";
    $numRemoved++;
    my $size = -s $file;
    $sizeRemoved += $size;
    if(!$$settings{debug}){
      unlink($file) or die "ERROR: could not unlink $file: $!";
    }
  }

  #print Dumper {size=>$sizeRemoved, num=>$numRemoved};
  return (size=>$sizeRemoved, num=>$numRemoved);
}

# https://stackoverflow.com/a/4493532
sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

sub usage{
  print "Cleans up a SneakerNet folder, usually to save some space

  Usage: $0 MiSeq_run_dir
  --debug  Just print what would be removed and do not remove
  \n";
  exit 0;
}

