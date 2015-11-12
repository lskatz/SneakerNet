#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;
use Email::Stuffer;
use List::MoreUtils qw/uniq/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

# All logging will go to a file, which will end up as $run/SneakerNet.txt.
# The link to the log will be emailed in the email plugin.
my $logdir=tempdir("$0_XXXXXX",TMPDIR=>1,CLEANUP=>1);
my $logfile="$logdir/logfile.txt";
open(my $logfileFh,'>',$logfile) or die "ERROR: could not open $logfile for writing: $!";
sub logmsg{
  my $msg="$0: @_\n";
  print STDERR $msg;
  print $logfileFh $msg;
}

exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug now test)) or die $!;
  die usage() if($$settings{help});

  if($$settings{test}){
    $$settings{now}=1;
    $$settings{inbox}=createTestDataset();
  } else {
    $$settings{inbox}||="/mnt/monolith0Data/dropbox/inbox";
  }

  # Find the directories
  my $inboxPath=$$settings{inbox};
  my $dirInfo=findReadsDir($inboxPath,$settings);
  # Process the directories
  for my $d(@$dirInfo){
    logmsg "Going to move $$d{dir}";
    waitForAnyChanges($d,$settings);
    #takeOwnership($d,$settings); # not sure if root really needs to own this first
    moveDir($d,$settings);
    giveToSequencermaster($d,$settings);

    # At this point, the log file should be put into this current directory.
    # Also, a SneakerNet directory should be created.
    # The file handle should be reopened.
    my $sneakernetDir="$$d{dir}/SneakerNet";
    mkdir $sneakernetDir;
    mkdir "$sneakernetDir/forEmail"; # anything in this directory will be emailed in the email plugin
    close $logfileFh;
    my $newLogfile="$sneakernetDir/SneakerNet.txt";
    system("cp -v $logfile $newLogfile"); # don't move, in case there is another run whose log needs to be copied
    die if $?;
    $logfile=$newLogfile;
    open($logfileFh,'>>',$logfile) or die "ERROR: could not open $logfile for writing: $!";
    symlink($logfile,"$$d{dir}/SneakerNet/forEmail/sneakernet.log");

    # Run all plugins as sequencermaster, using ssh to act as sequencermaster
    my @exe=glob("$FindBin::RealBin/SneakerNet.plugins/*");

    # Put the email script in last and remove files that are
    # not executable and not files.
    @exe=grep {!/emailWhoever.pl/} @exe;
    push(@exe,"$FindBin::RealBin/SneakerNet.plugins/emailWhoever.pl");
    @exe=map{ 
      $_ if(-f $_ && -x $_);
    } @exe;
    @exe=grep {/./} @exe;

    # Execute all plugins
    for(my $i=0;$i<@exe;$i++){
      command("ssh sequencermaster\@localhost $exe[$i] $$d{dir} --numcpus $$settings{numcpus}");
    }
  }

  return 0;
}

sub createTestDataset{
  my $inbox=tempdir("SneakerNetXXXXXX",TMPDIR=>1,CLEANUP=>1);
  my $rundir="$inbox/test-15-001";
  mkdir($rundir);

  # create fastq files
  for my $sample(qw(A B)){
    my $forward="$rundir/${sample}_1.fastq";
    my $reverse="$rundir/${sample}_2.fastq";
    logmsg "Creating $forward, $reverse";
    my $read="A" x 150; 
    my $qual="I" x 150;
    open(FWD,">",$forward) or die "ERROR: could not make temp file $forward: $!";
    open(REV,">",$reverse) or die "ERROR: could not make temp file $reverse: $!";
    for my $i(1..4e5){
      # read1
      print FWD "\@".$i."/1\n$read\n+\n$qual\n";
      # read2
      print REV "\@".$i."/2\n$read\n+\n$qual\n";
    }
    close FWD;
    close REV;
    command("gzip $forward $reverse");
  }

  # create spreadsheet
  my $samplesheet="$rundir/SampleSheet.csv";
  open(SAMPLE,">",$samplesheet) or die "ERROR: cannot write to $samplesheet: $!";
  print SAMPLE "[Header]\nIEMFileVersion,4\nInvestigator Name,LSK (gzu2)\nExperiment Name,test\n";
  print SAMPLE "Date,8/14/2015\nWorkflow,GenerateFASTQ\nApplication,FASTQ Only\nAssay,Nextera XT\nDescription,Listeria GMI\nChemistry,Amplicon\n\n";
  print SAMPLE "[Reads]\n150\n150\n]n";
  print SAMPLE "[Settings]\nReverseComplement,0\nAdapter,CTGTCTCTTATACACATCT\n\n";
  print SAMPLE "[Data]\nSample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description\n";
  print SAMPLE "A,,test,A01,N701,TAAGGCGA,S517,GCGTAAGA,,Species=Listeria_monocytogenes;ExpectedGenomeSize=0.4;Route=CalcEngine;Route=NCBI\n";
  print SAMPLE "B,,test,B01,N702,CGTACTAG,S517,GCGTAAGA,,Species=Listeria_monocytogenes;ExpectedGenomeSize=0.4;Route=NCBI\n";
  close SAMPLE;

  # make zero byte files
  for my $i(qw(config.xml SampleSheet.csv QC/CompletedJobInfo.xml QC/InterOp/CorrectedIntMetricsOut.bin QC/runParameters.xml QC/GenerateFASTQRunStatistics.xml QC/RunInfo.xml)){
    my $zerobyte="$rundir/$i";
    mkdir dirname($zerobyte);
    logmsg "Making zero-byte file $zerobyte";
    open(FILE,">>", $zerobyte) or die "ERROR: could not create zero-byte file $zerobyte: $!";
    print FILE "blah\n";
    close FILE;
  }

  return $inbox;
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
     $waitSeconds=0 if($$settings{now});
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
  my $destinationDir="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
  if(-e $destinationDir){
    die "ERROR: destination directory already exists!\n  $destinationDir";
  }
  command("mv --no-clobber -v $$info{dir} $destinationDir");

  $$info{source_dir}=$$info{dir};
  $$info{subdir}=$subdir;
  $$info{dir}="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
}

sub giveToSequencermaster{
  my($info,$settings)=@_;
  command("chown -Rv sequencermaster.sequencermaster $$info{dir}");
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
  my $stdout=`$command 2>&1`;
  my $exit_code=$?;
  print STDERR $stdout;
  print $logfileFh $stdout;
  die "ERROR running command\n  $command" if $exit_code;
}

sub usage{
  "Find all reads directories under the inbox, puts them into the right
  place on Monolith0 and gives ownership to sequencermaster.
  All executable scripts under SneakerNet.plugins will also be run.

  Usage: $0 [-i inboxDir/]
  -i dir  # choose a different 'inbox' to look at
  --test  # Create a test directory 
  --debug # Show debugging information
  --now   # Get this show on the road!!
  "
}
