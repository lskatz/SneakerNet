#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help inbox=s debug));
  die usage() if($$settings{help});
  my $inboxPath=$$settings{inbox}||"/media/data/dropbox/inbox";

  # Find the directories
  my $dirInfo=findReadsDir($inboxPath,$settings);
  # Process the directories
  for my $d(@$dirInfo){
    takeOwnership($d,$settings);
    moveDir($d,$settings);
    addReadMetrics($d,$settings);
    giveToSequencermaster($d,$settings);
  }

  return 0;
}

sub findReadsDir{
  my($inbox,$settings)=@_;

  # all subdirectories under the inbox
  my @dir=grep({-d $_} glob("$inbox/*"));

  my @readsDir; # holds the true reads directories
  for my $d(@dir){
    my $info=parseReadsDir($d,$settings);
    logmsg Dumper $info if($$settings{debug});
    push(@readsDir,$info) if($$info{is_good});
  } 

  return \@readsDir;
}

# Figure out if this really is a reads directory
sub parseReadsDir{
  my($dir,$settings)=@_;

  my %dirInfo=(dir=>$dir,is_good=>1,why_not=>"");

  my $b=basename $dir;
  ($dirInfo{machine},$dirInfo{year},$dirInfo{run},$dirInfo{comment})=split(/\-/,$b);

  # If the run name isn't even there, then it's not a run directory
  if(!defined($dirInfo{run})){
    $dirInfo{why_not}.="Run is not defined for $dir\n";
    $dirInfo{is_good}=0;
  }

  # See if there are actually reads in the directory
  if(!glob("$dir/*.fastq.gz")){
    $dirInfo{why_not}.= "Could not find fastq.gz files in $dir\n";
    $dirInfo{is_good}=0;
  }

  # See if the misc. files are in there too
  for(qw(config.xml SampleSheet.csv QC/CompletedJobInfo.xml QC/InterOp QC/runParameters.xml QC/GenerateFASTQRunStatistics.xml QC/RunInfo.xml)){
    if(!-e "$dir/$_"){
      $dirInfo{why_not}.="Could not find $dir/$_\n";
      $dirInfo{is_good}=0;
    }
  }

  return \%dirInfo;
}

sub takeOwnership{
  my($info,$settings)=@_;
  my $user=$ENV{USER};
  command("chown -R $user.$user $$info{dir}");
  die if $?;
}

sub moveDir{
  my($info,$settings)=@_;

  my $subdir=join("-",$$info{machine},$$info{year},$$info{run});
  command("mv --no-clobber -v $$info{dir} /media/data/RawSequenceData/$$info{machine}/$subdir");
  $$info{source_dir}=$$info{dir};
  $$info{dir}="/media/data/RawSequenceData/$$info{machine}/$subdir";
}

sub addReadMetrics{
  my($info,$settings)=@_;
  command("run_assembly_readMetrics.pl --fast $$info{dir}/*.fastq.gz | sort -k3,3n > $$info{dir}/readMetrics.txt");
}
 
sub giveToSequencermaster{
  my($info,$settings)=@_;
  command("chown -R sequencermaster.sequencermaster $$info{dir}");
}


sub command{
  my($command,$settings)=@_;
  logmsg "COMMAND\n  $command" if($$settings{debug});
  system($command);
  die "ERROR running command\n  $command" if $?;
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 [-i inboxDir/]
  -i dir  # choose a different 'inbox' to look at
  --debug # Show debugging information
  "
}
