#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy mv cp/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;
use List::MoreUtils qw/uniq/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig/;

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
  GetOptions($settings,qw(help numcpus=i inbox=s debug now test email! preserve)) or die $!;
  die usage() if($$settings{help});
  $$settings{numcpus}||=1;
  $$settings{email}//=1;
  $$settings{preserve}//=0;

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

    moveDir($d,$settings);

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
    link($logfile,"$$d{dir}/SneakerNet/forEmail/sneakernet.log"); # hard link to help resolve relative path issues

    # also add in the Sample Sheet for email for later
    link("$$d{dir}/SampleSheet.csv","$sneakernetDir/forEmail/SampleSheet.csv");

    # Give the rest to sequencermaster, now that it has all been moved over
    # Deprecated: All sequences are now copied over by sequencermaster and are owned by sequencermaster.
    # system("chown -R sequencermaster.sequencermaster $$d{dir}/SneakerNet");

    #my @exe=flatten($$settings{'plugins.conf'});
    my @exe=@{ $$settings{'default.plugins'} };
    for my $exe(@exe){
      if(!$$settings{email} && $exe=~/emailWhoever.pl/){
        next;
      }

      command("$FindBin::RealBin/../SneakerNet.plugins/$exe $$d{dir} --numcpus $$settings{numcpus} 2>&1");
    }
    
    # Add permissions for the sequencermaster group
    command("chmod -R g+wr $$d{dir}");
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

  my %dirInfo=(dir=>$dir,is_good=>1,why_not=>"", runType=>"");

  my $b=basename $dir;
  ($dirInfo{machine},$dirInfo{year},$dirInfo{run},$dirInfo{comment})=split(/\-/,$b);

  # If the run name isn't even there, then it's not a run directory
  if(!defined($dirInfo{run})){
    $dirInfo{why_not}.="Run name is not defined for $dir. Run name syntax should be Machine-year-runNumber-comment.\n";
    $dirInfo{is_good}=0;
    return \%dirInfo;
  }

  # Test for Illumina at the same time as seeing if all the files are in there
  if(!$dirInfo{is_good} || !$dirInfo{runType}){

    my $foundAllFiles=1;

    # See if there are actually reads in the directory
    if(!glob("$dir/*.fastq.gz")){
      $dirInfo{why_not}.= "[Illumina] Could not find fastq.gz files in $dir\n";
      $foundAllFiles=0;
    }

    # How do we tell it is a miniseq run?  My best guess
    # is if we see "SampleSheetUsed.csv" instead of
    # "SampleSheet.csv."
    if(-e "$dir/SampleSheetUsed.csv"){
      logmsg "Detected $dir/SampleSheetUsed.csv: it could be a miniseq run.";
      # cp the sample sheet to SampleSheet.csv to make it compatible.
      cp("$dir/SampleSheetUsed.csv","$dir/SampleSheet.csv");
      cp("$dir/QC/RunParameters.xml","$dir/QC/runParameters.xml");

      # edit the sample sheet to remove the run
      removeRunNumberFromSamples("$dir/SampleSheet.csv", $settings);
      
      # Make empty files for compatibility
      for("$dir/config.xml"){
        open(EMPTYFILE,">>", $_) or die "ERROR: could not make an empty file $_: $!";
        close EMPTYFILE;
      }
    }

    # See if the misc. files are in there too
    for(qw(config.xml SampleSheet.csv QC/CompletedJobInfo.xml QC/InterOp QC/runParameters.xml QC/GenerateFASTQRunStatistics.xml QC/RunInfo.xml)){
      if(!-e "$dir/$_"){
        $dirInfo{why_not}.="[Illumina] Could not find $dir/$_\n";
        $foundAllFiles=0;
      }
    }
    
    $dirInfo{runType}="Illumina" if($foundAllFiles);
  }

  # Test for Ion Torrent at the same time as seeing if all the files are in there
  if(!$dirInfo{is_good} || !$dirInfo{runType}){

    my $foundAllFiles=1;
    
    # See if there are reads in the directory
    my @fastq=(glob("$dir/plugin_out/downloads/*.fastq"),glob("$dir/plugin_out/downloads/*.fastq.gz"));
    if(!@fastq){
      $dirInfo{why_not}.= "[IonTorrent] Could not find fastq[.gz] files in $dir/plugin_out/downloads\n";
      $foundAllFiles=0;
    }
    
    $dirInfo{runType}="IonTorrent" if($foundAllFiles);
  }

  if(!$dirInfo{runType}){
    system("mv -v $dir $$settings{inbox}/rejected");
    logmsg "WARNING: tried to move $dir into $$settings{inbox}/reject but was not able to" if($?);
    logmsg "ERROR: could not determine the run type of $dir (e.g., Illumina or IonTorrent). Additional info to complete the run for any particular chemistry:\n$dirInfo{why_not}";
  }

  #die Dumper \%dirInfo;


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

sub moveDir{
  my($info,$settings)=@_;

  $$info{comment}||="";
  my $subdir=join("-",$$info{machine},$$info{year},$$info{run},$$info{comment});
  $subdir=~s/\-$//; # remove final dash in case the comment wasn't there
  my $destinationDir="$$settings{REPOSITORY_DIRECTORY}/$$info{machine}/$subdir";
  if(-e $destinationDir){
    die "ERROR: destination directory already exists!\n  $destinationDir";
  }

  #die Dumper $info;
  # Copy and then delete, so that permissions are retained for sequencermaster
  command("mkdir -p $$settings{REPOSITORY_DIRECTORY}/$$info{machine}");
  command("cp --no-clobber -vr $$info{dir} $destinationDir");
  command("rm -vfr $$info{dir}") if(!$$settings{preserve});

  # Add permissions for the sequencermaster group
  command("chmod -R g+wr $destinationDir");

  # Update some attributes about this run
  $$info{source_dir}=$$info{dir};
  $$info{subdir}=$subdir;
  $$info{dir}=$destinationDir;

  # Now that the source directory has been 'moved' to a 
  # new location, in the event of an error, it should
  # be moved back to the rejects folder.
  $SIG{__DIE__} = sub{
    my $rejectFolder="$$settings{inbox}/rejected";
    mkdir $rejectFolder;
    chmod(oct("2775"),$rejectFolder); # drwxrwsr-x

    # Don't use command() because it has a potential die
    # statement in there and could cause an infinite loop.
    eval{
      system("cp -rv $destinationDir $rejectFolder/");
      system("rm -rf $destinationDir");
      logmsg "ERROR: Moved the error folder from $destinationDir to $rejectFolder";
    };
    if($@){
      logmsg "ERROR: there was a problem copying $destinationDir to $rejectFolder/.";
    }

    die @_;
  };
}

# Edit a sample sheet in-place to remove a run identifier
# from the sample names. For some reason the Miniseq
# appends a four digit number, e.g. "-6006" to the end
# of each sample name.
sub removeRunNumberFromSamples{
  my($samplesheet,$settings)=@_;

  my $newSamplesheetString="";
  open(SAMPLESHEET,"<", $samplesheet) or die "ERROR: could not read $samplesheet: $!";
  my $reachedSamples=0;
  my $runid="";
  while(<SAMPLESHEET>){
    # Make a note of the run ID when I see it
    if(/Local Run Manager Analysis Id,\s*(\d+)/){
      $runid=$1;
    }

    if(!$reachedSamples){
      $newSamplesheetString.=$_;
      if(/Sample_ID,/){
        $reachedSamples=1;
      }
    }
    # Read the samples and remove the run ID
    else {
      my($samplename,@therest)=split(/,/,$_);
      $samplename=~s/\-$runid$//;
      $newSamplesheetString.=join(",",$samplename,@therest);
    }
  }
  close SAMPLESHEET;

  # Now rewrite the sample sheet
  open(SAMPLESHEET,">", $samplesheet) or die "ERROR: could not write to $samplesheet: $!";
  print SAMPLESHEET $newSamplesheetString;
  close SAMPLESHEET;

  return 1;
}

################
# Utility subs #
################

sub readConfigOld{
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
  logmsg "Executing command\n    $command";
  my $stdout=SneakerNet::command($command,$settings);

  print $logfileFh $stdout;

  return $stdout;
}

sub flatten {
  map { ref $_ ? flatten(@{$_}) : $_ } @_;
}

sub usage{
  "Find all reads directories under the inbox, puts them into the right
  place on Monolith0 and gives ownership to sequencermaster.
  All executable scripts under SneakerNet.plugins will also be run.

  Usage: $0 [-i inboxDir/]
  --inbox  dir  # choose a different 'inbox' to look at
  --test        # Create a test directory. Implies --now.
  --noemail     # Do not send an email at the end.
  --preserve    # Do not delete the source run (Default: cp and then rm -r)
  --debug       # Show debugging information
  --now         # Do not check whether the directory contents are still being modified.
  --numcpus  1
  "
}
