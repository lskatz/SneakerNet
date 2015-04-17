#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use FindBin;
use MIME::Lite;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug force));
  die usage() if($$settings{help});
  $$settings{inbox}||="/mnt/monolith0Data/dropbox/inbox";

  # Find the directories
  my $inboxPath=$$settings{inbox};
  my $dirInfo=findReadsDir($inboxPath,$settings);
  # Process the directories
  for my $d(@$dirInfo){
    logmsg "Going to move $$d{dir}";
    waitForAnyChanges($d,$settings);
    takeOwnership($d,$settings);
    moveDir($d,$settings);
    addReadMetrics($d,$settings);
    giveToSequencermaster($d,$settings);
    emailWhoever($d,$settings);
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
    if($$info{is_good}){
      push(@readsDir,$info) if($$info{is_good});
    } elsif($$info{why_not}){
      logmsg "Found $d but $$info{why_not}";
    }
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
  $waitSeconds=0 if($$settings{force});
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

  $$info{comment}||="";
  my $subdir=join("-",$$info{machine},$$info{year},$$info{run},$$info{comment});
  $subdir=~s/\-$//; # remove final dash in case the comment wasn't there
  command("mv --no-clobber -v $$info{dir} /mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir");

  $$info{source_dir}=$$info{dir};
  $$info{subdir}=$subdir;
  $$info{dir}="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
}

sub addReadMetrics{
  my($info,$settings)=@_;
  logmsg "Running fast read metrics";
  command("run_assembly_readMetrics.pl --fast $$info{dir}/*.fastq.gz | sort -k3,3n > $$info{dir}/readMetrics.txt.tmp");

  # edit read metrics to include genome sizes
  my $newReadMetrics;
  open(READMETRICS,"$$info{dir}/readMetrics.txt.tmp") or die "ERROR: could not open $$info{dir}/readMetrics.txt.tmp because $!";
  open(READMETRICSFINAL,">","$$info{dir}/readMetrics.txt") or die "ERROR: could not open $$info{dir}/readMetrics.txt for writing: $!";

  # get the header and also put it into the final output file
  my $header=<READMETRICS>;
  print READMETRICSFINAL $header;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<READMETRICS>){
    chomp;
    # read in each line into the appropriate header
    my %h;
    @h{@header}=split(/\t/);

    # find the genome size based on the filename
    my $coverage=calculateCoverage(\%h,$settings);
    $h{coverage}=$coverage;

    for(@header){
      print READMETRICSFINAL "$h{$_}\t";
    }
    print READMETRICSFINAL "\n";
  }
  close READMETRICSFINAL;
  close READMETRICS;
}

# Use a hash of a line from a readmetrics output 
# to determine genome coverage.
sub calculateCoverage{
  my($h,$settings)=@_;
  my $coverage=$$h{coverage} || '.'; # default value in case one isn't found
  my $is_recalculated=0;             # know whether the coverage was recalculated

  # See if we can recalculate the coverage based on the filename
  my $file=basename($$h{File});
  for my $info(@{ $$settings{genomeSizes} }){
    my($regex,$size,$organism)=@$info;
    next if($file!~/$regex/);

    $coverage=$$h{totalBases}/$size;
    $coverage=sprintf("%0.2f",$coverage); # round it
    $is_recalculated=1; # the coverage was recalculated

    logmsg "Decided that $$h{File} is $organism with expected genome size $size. New coverage: $coverage";
  }

  # Report!
  if(!$is_recalculated){
    logmsg "Warning: could not understand what organism $$h{File} belongs to; coverage was not recalculated. Reported coverage: $coverage";
  }
  return $coverage;
}
 
sub giveToSequencermaster{
  my($info,$settings)=@_;
  command("chown -Rv sequencermaster.sequencermaster $$info{dir}");
}

sub emailWhoever{
  my($info,$settings)=@_;
  
  my $subdir=$$info{subdir};
  my $readMetrics=$$info{dir}."/readMetrics.txt";

  my $from="sequencermaster\@monolith0.edlb.cdc.gov";
  my $to="gzu2\@cdc.gov";
  my $uuencodeAttach=`uuencode $subdir.qc.txt < $readMetrics`;

  # This worked for me:
  # (echo -e "Subject: M947-15-012 QC\nFrom: root@monolith0.edlb.cdc.gov"; uuencode blah.txt < /mnt/monolith0Data/RawSequenceData/M947/M947-15-012/readMetrics.txt) | sendmail gzu2@cdc.gov

  logmsg "Emailing to $to\n  readMetrics file: $readMetrics";

  my $fullMessage="Subject: $subdir QC\nFrom: $from\nTo: $to\nPlease open the following attachment in Excel for read metrics for run $subdir.\n$uuencodeAttach";

  my $exit_code=system("echo -e \"$fullMessage\" | sendmail $to");
  return !$exit_code;
}

################
# Utility subs #
################
sub readConfig{
  my @file=glob("$FindBin::RealBin/config/*");
  my $settings={};
  for(@file){
    open(CONFIGFILE,$_) or die "ERROR: could not open config file $_: $!";
    my $key=basename $_;
    while(<CONFIGFILE>){
      s/^\s+|\s+$//g; # trim
      next if(/^$/);
      next if(/^#/);
      my $configLine=[split(/\t/,$_)];
      push(@{ $$settings{$key} },$configLine);
    }
    close CONFIGFILE;
  }
  return $settings;
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
  --force # Get this show on the road!!
  "
}
