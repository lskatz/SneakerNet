#!/usr/bin/env perl
# Transfers files to a remote destination and QCs them beforehand.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use List::Util qw/sum/;

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig samplesheetInfo command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  transferFilesToRemoteComputers($dir,$settings);

  return 0;
}

sub transferFilesToRemoteComputers{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  # Which files should be skipped according to Q/C?
  my $tmp=identifyBadRuns($dir,$sampleInfo,$settings);
  my $toSkip=$$tmp{toSkip};
  my $whatFailed=$$tmp{whatFailed};

  # Make a file detailing what passed or failed.
  my @header=qw(File coverage quality failed);
  my $passfail="$dir/SneakerNet/forEmail/passfail.tsv";
  open(PASSFAIL,">",$passfail) or die "ERROR: could not open $passfail for writing: $!";
  print PASSFAIL join("\t",@header)."\n";
  
  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases
    my $taxon=$$s{species} || 'NOT LISTED';
    logmsg "The taxon of $sampleName is $taxon";
    if(grep {/calcengine/i} @{ $$s{route} }){
      for(@{ $$s{fastq} }){
        # Write out the status
        # The key of each filename is its basename
        my $f=basename($_);
        # Set a default hash here
        $$whatFailed{$f}//={coverage=>0,quality=>0};

        # The boolean $failed is calculated based on whether any of the values are true
        # and should be equal to what is in %$toSkip.
        my $failed=int(!!sum(values(%{ $$whatFailed{$f} })));
        $$whatFailed{$f}{coverage}//=0;
        $$whatFailed{$f}{quality}//=0;
        print PASSFAIL join("\t",$f,$$whatFailed{$f}{coverage},$$whatFailed{$f}{quality},$failed)."\n";

        next if($$toSkip{basename($_)});
        my $subfolder=$$s{taxonRules}{dest_subfolder} || "SneakerNet";
        $filesToTransfer{$subfolder}.=$_." ";
      }
      logmsg "One route for sample $sampleName is the Calculation Engine";
    } else {
      logmsg "Note: The route for $sampleName was not listed in the sample sheet.";
    }
  }
  print PASSFAIL "# 1 (one) indicates a failure in a particular category.\n";
  print PASSFAIL "# 1 in the failed category means that the sample failed in any category.\n";
  close PASSFAIL;

  #die "ERROR: no files to transfer" if (!$filesToTransfer);
  logmsg "WARNING: no files will be transferred" if(!keys(%filesToTransfer));

  # Make the transfers based on taxon.
  while(my($subfolder,$fileString)=each(%filesToTransfer)){

    logmsg "Transferring to $subfolder:\n  $fileString";
    if(!$$settings{debug}){
      command("rsync --update -av --no-g $fileString $$settings{transfer_destination_string}/$subfolder/");
    }
  }
}

# Which runs should be skipped based on bad quality
# or bad coverage, or whatever?
# Returns filenames.
sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %toSkip=();      # boolean fail
  my %whatFailed=();  # reasons why it failed

  open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
  my @header=split(/\t/,<READMETRICS>); chomp(@header);
  while(<READMETRICS>){
    chomp;
    my %F;
    @F{@header}=split(/\t/,$_);
    my $samplename=basename($F{File},'.fastq.gz');
    $samplename=~s/_S\d+_.*//; # figure out the sample name before the the _S1_ pattern

    # Skip anything that says undetermined.
    if($samplename=~/^Undetermined/){
      $toSkip{basename($F{File})}=1;
      next;
    }

    #die Dumper [$samplename,$$sampleInfo{$samplename},\%F];

    # Get the name of all files linked to this file through the sample.
    my @file=@{$$sampleInfo{$samplename}{fastq}};
    
    # Compare coverage of one read against half of the
    # threshold coverage because of PE reads.
    if($F{coverage} ne '.'){ # dot means coverage is unknown.
      if($F{coverage} < $$sampleInfo{$samplename}{taxonRules}{coverage}/2){
        for (@file){
          my $f=basename($_); # avoid the directory name
          $toSkip{$f}=1;
          $whatFailed{$f}{coverage}=1;
        }
        logmsg "Low coverage in $F{File}";
      }
    }

    if($F{avgQuality} < $$sampleInfo{$samplename}{taxonRules}{quality}){
      for (@file){
        my $f=basename($_); # avoid the directory name
        $toSkip{$f}=1;
        $whatFailed{$f}{avgQuality}=1;
      }
      logmsg "low quality in $F{File}\n  Skipping @file";
    }

  }

  return {toSkip=>\%toSkip,whatFailed=>\%whatFailed};
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 MiSeq_run_dir
  "
}

