#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use Email::Stuffer;
use List::MoreUtils qw/uniq/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug force test)) or die $!;
  die usage() if($$settings{help});

  if($$settings{test}){
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
    takeOwnership($d,$settings);
    moveDir($d,$settings);
    addReadMetrics($d,$settings);
    giveToSequencermaster($d,$settings);
    transferFilesToRemoteComputers($d,$settings);
    emailWhoever($d,$settings);
  }

  return 0;
}

sub createTestDataset{
  my $inbox=File::Temp->tempdir("SneakerNetXXXXXX",TMPDIR=>1,CLEANUP=>1);
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
  my $destinationDir="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
  if(-e $destinationDir){
    die "ERROR: destination directory already exists!\n  $destinationDir";
  }
  command("mv --no-clobber -v $$info{dir} $destinationDir");

  $$info{source_dir}=$$info{dir};
  $$info{subdir}=$subdir;
  $$info{dir}="/mnt/monolith0Data/RawSequenceData/$$info{machine}/$subdir";
}

sub addReadMetrics{
  my($info,$settings)=@_;
  logmsg "Running fast read metrics";
  command("run_assembly_readMetrics.pl --fast $$info{dir}/*.fastq.gz | sort -k3,3n > $$info{dir}/readMetrics.tsv.tmp");

  # edit read metrics to include genome sizes
  my $newReadMetrics;
  open(READMETRICS,"$$info{dir}/readMetrics.tsv.tmp") or die "ERROR: could not open $$info{dir}/readMetrics.tsv.tmp because $!";
  open(READMETRICSFINAL,">","$$info{dir}/readMetrics.tsv") or die "ERROR: could not open $$info{dir}/readMetrics.tsv for writing: $!";

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

  # Clean up by removing the temporary file
  unlink("$$info{dir}/readMetrics.tsv.tmp");
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

sub transferFilesToRemoteComputers{
  my($info,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$$info{dir}/SampleSheet.csv",$settings);

  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    my $taxon=(keys(%{ $$s{species} }))[0];
    if($$s{route}{calcengine}){
      $filesToTransfer{$taxon}.=join(" ",glob("$$info{dir}/$sampleName*.fastq.gz"))." ";
    }
  }

  #die "ERROR: no files to transfer" if (!$filesToTransfer);
  logmsg "WARNING: no files will be transferred" if(!keys(%filesToTransfer));

  # Make the transfers based on taxon.
  # TODO consider putting this taxon logic into a config file.
  while(my($taxon,$fileString)=each(%filesToTransfer)){

    # Which folder under /scicomp/groups/OID/NCEZID/DFWED/EDLB/share/out/Calculation_Engine
    # is appropriate?  SneakerNet if nothing else is found.
    my $subfolder="SneakerNet";
    if($taxon =~ /Listeria/i){
      $subfolder="LMO";
    } elsif ($taxon =~ /Salmonella/i){
      $subfolder="Salm";
    } elsif ($taxon =~ /Campy|Arcobacter|Helicobacter/i){
      $subfolder="Campy";
    } elsif ($taxon =~ /^E\.$|STEC|Escherichia|Shigella/i){
      $subfolder="STEC";
    } elsif ($taxon =~ /Vibrio|cholerae|cholera/i){
      $subfolder="Vibrio";
    } else {
      logmsg "WARNING: cannot figure out the correct subfolder for taxon $taxon. The following files will be sent to $subfolder instead.";
    }
    logmsg "Transferring to $subfolder:\n  $fileString";
    command("rsync --update -av $fileString gzu2\@aspen.biotech.cdc.gov:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/out/Calculation_Engine/$subfolder/");
  }
}

sub samplesheetInfo{
  my($samplesheet,$settings)=@_;

  my $section="";
  my @header=();
  my %sample;
  open(SAMPLE,$samplesheet) or die "ERROR: could not open sample spreadsheet $samplesheet: $!";
  while(<SAMPLE>){
    chomp;
    if(/^\[(\w+)\]/){
      $section=lc($1);
      my $header=<SAMPLE>;
      chomp($header);
      @header=split(/,/,lc($header));
      next;
    }
    if($section eq "data"){
      my %F;
      @F{@header}=split(/,/,$_);
      for my $keyvalue(split(/;/,lc($F{description}))){
        my($key,$value)=split(/=/,$keyvalue);
        $key=~s/^\s+|\s+$//g;      #whitespace trim
        $value=~s/^\s+|\s+$//g;    #whitespace trim
        $F{$key}={} if(!$F{$key});
        $F{$key}{$value}++;
      }
      
      $sample{$F{sample_id}}=\%F;
    }
  }
  return \%sample;
}

sub emailWhoever{
  my($info,$settings)=@_;

  my $subdir=$$info{subdir};
  my $readMetrics=$$info{dir}."/readMetrics.tsv";

  # Figure out who we are mailing
  my @to=("gzu2\@cdc.gov","wwm8\@cdc.gov","pfge\@cdc.gov","wvt2\@cdc.gov","fid4\@cdc.gov");
  my $pocLine=`grep -m 1 'Investigator' $$info{dir}/SampleSheet.csv`;
  if($pocLine=~/\((.+)\)/){
    my $cdcids=$1;
    $cdcids=~s/\s+//g;
    my @cdcids=map {"$_\@cdc.gov"} split(/,/,$cdcids);
    push(@to,@cdcids);
  } else {
    logmsg "WARNING: could not parse the investigator line so that I could find CDC IDs";
  }

  @to=uniq(@to);

  # Send one email per recipient.
  for my $to(@to){
    logmsg "To: $to";
    my $from="sequencermaster\@monolith0.edlb.cdc.gov";
    my $subject="$subdir QC";
    my $body ="Please open the following attachment in Excel for read metrics for run $subdir.\n";
       $body.="\n  This message was brought to you by SneakerNet!";

    my $was_sent=Email::Stuffer->from($from)
                               ->subject($subject)
                               ->to($to)
                               ->text_body($body)
                               ->attach_file($readMetrics)
                               ->send;

    if(!$was_sent){
      logmsg "Warning: Email was not sent to $to!";
    }
  }

  return \@to;
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
  --test  # Create a test directory 
  --debug # Show debugging information
  --force # Get this show on the road!!
  "
}
