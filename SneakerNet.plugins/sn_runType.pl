#!/usr/bin/env perl
# Figure out what kind of run this is

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use File::Spec;
use Text::Fuzzy;

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig logmsg samplesheetInfo /;
use Email::Stuffer;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help force tempdir=s debug numcpus=i email:s )) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  #logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  my $dirInfo = parseReadsDir($dir,$settings);

  my $runStatus="";
  my $exitStatus=0;

  if($$dirInfo{runType}){
    print $$dirInfo{runType}."\n";
    $runStatus.="Run type for $$dirInfo{run_name} is $$dirInfo{runType}, and it is ready to run!\n";
    
  } else {
    $runStatus.="ERROR: could not determine the run type of $dir. Additional info to complete the run for any particular chemistry:\n$$dirInfo{why_not}\n";
    $exitStatus=3;
  }

  if($$settings{debug}){
    logmsg $runStatus;
  }

  # Send email using snok.txt?
  if(defined(my $to = $$settings{email})){

    # append any snok.txt emails
    # Read the run's snok.txt for any emails
    if(-e "$dir/snok.txt"){
      my @email;
      my $snokCfg = new Config::Simple();
      eval{
        $snokCfg->read("$dir/snok.txt");
        @email = $snokCfg->param("emails");
        logmsg "Found emails in snok.txt: @email";
      };
      if($@){
        logmsg "WARNING: could not read snok.txt for any emails, but the file exists!";
      }
      $to.=",".join(",",@email);
    }
      
    my $from=$$settings{from} || die "ERROR: need to set 'from' in the settings.conf file!";
    my $subject="Initial SneakerNet status for $$dirInfo{run_name}";
    my $body="SneakerNet took a first glance at $$dirInfo{run_name} and is reporting...\n$runStatus";
    my $email=Email::Stuffer->from($from)
                            ->subject($subject)
                            ->to($to)
                            ->text_body($body);
    $email->send;
    logmsg "Email sent!";
  }
  logmsg "Exit status: $exitStatus";

  return $exitStatus;
}

# Figure out if this really is a reads directory
sub parseReadsDir{
  my($dir,$settings)=@_;

  my %dirInfo=(dir=>$dir,is_good=>1,why_not=>"", runType=>"");

  # The directory name, regardless 
  my $b=basename(File::Spec->rel2abs($dir));
  ($dirInfo{machine},$dirInfo{year},$dirInfo{run},$dirInfo{comment})=split(/\-/,$b);
  $dirInfo{run_name}=$b;

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

  # Now that we know the run type, see if the sample sheet
  # makes sense.
  if($dirInfo{runType} eq "Illumina"){
    my @fastq = glob("$dir/*.fastq.gz");
    my $samples = samplesheetInfo("$dir/SampleSheet.csv",$settings);
    while(my($sample,$sampleHash) = each(%$samples)){
      next if(ref($sampleHash) ne 'HASH');
      my $tf = Text::Fuzzy->new($sample);

      # Even if there are fastq files for this sample, double
      # check that they exist.
      for my $candidateFastq(@{$$sampleHash{fastq}}){

        if(-e $candidateFastq){
          next;
        }

        my $nearest = $tf->nearestv(\@fastq);
        $dirInfo{why_not}.="For sample $sample, it looks like the fastq file $candidateFastq does not exist.  Did you mean $nearest?\n";
        $dirInfo{is_good}=0;
        $dirInfo{runType}="";
      }

      # If there are zero fastq files for this sample, make
      # suggestions
      if(!@{$$sampleHash{fastq}}){
        my $nearest = $tf->nearestv(\@fastq);
        $dirInfo{why_not}.="Could not find exact matches for fastq files for sample $sample. For this sample, the closest file match is $nearest\n";
        $dirInfo{is_good}=0;
        $dirInfo{runType}="";
      }
    }
        
  }

  return \%dirInfo;
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


sub usage{
  "Print the type of run directory.\nExit code 3 if the run is invalid.

  Usage: $0 [options] -- dir/
  --debug   Print additional information about the run to stderr
  --email   If supplied, even without a value, snok.txt will be
            used for report recipients.  If a value is given,
            the report will also be sent to that email.
  "
}

