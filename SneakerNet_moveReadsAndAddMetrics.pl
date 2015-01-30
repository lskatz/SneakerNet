#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help inbox=s debug));
  die usage() if($$settings{help});
  my $inboxPath=$$settings{inbox}||"/mnt/monolith0Data/dropbox/inbox";

  # Find the directories
  my $dirInfo=findReadsDir($inboxPath,$settings);
  # Process the directories
  for my $d(@$dirInfo){
    logmsg "Going to move $$d{dir}";
    waitForAnyChanges($d,$settings);
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

sub waitForAnyChanges{
  my($info,$settings)=@_;

  # Find a unique md5sum for this timepoint.
  # If it changes, then the directory has changed
  my $md5sumCommand="find $$info{dir} -type f -exec ls -l {} \\; | sort | md5sum -";
  my $md5sum=`$md5sumCommand`;
  die "ERROR with md5sum command: $!\n  $md5sum" if $?;
  
  my $waitSeconds=120; # two minutes just in case
  logmsg "I will wait $waitSeconds seconds to see if any changes are in progress in the directory.";
  for(my $i=1;$i<=$waitSeconds;$i++){
    sleep 1;
    logmsg "$i seconds..." if($i % 30 == 0);
    my $newMd5=`$md5sumCommand`;
    die "ERROR with md5sum command: $!\n  $md5sum" if $?;
    if($newMd5 ne $md5sum){
      logmsg "WARNING: directory has changed!  I will wait for another $waitSeconds seconds before I test again";
      $md5sum=$newMd5;
      $i=1;
    }
  }
  logmsg "No changes detected. Onward!";
  return 1;
}

sub takeOwnership{
  my($info,$settings)=@_;
  my $user=$ENV{USER};
  command("chown -Rv $user.$user $$info{dir}");
  die if $?;
}

sub moveDir{
  my($info,$settings)=@_;

  my $subdir=join("-",$$info{machine},$$info{year},$$info{run});
  command("mv --no-clobber -v $$info{dir} /mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir");
  $$info{source_dir}=$$info{dir};
  $$info{dir}="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
}

sub addReadMetrics{
  my($info,$settings)=@_;
  logmsg "Running fast read metrics";
  command("run_assembly_readMetrics.pl --fast $$info{dir}/*.fastq.gz | sort -k3,3n > $$info{dir}/readMetrics.txt");
}
 
sub giveToSequencermaster{
  my($info,$settings)=@_;
  command("chown -Rv sequencermaster.sequencermaster $$info{dir}");
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
