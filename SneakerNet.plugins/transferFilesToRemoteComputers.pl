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
use POSIX qw/strftime/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg passfail/;

our $VERSION = "1.7.1";
our $CITATION= "Transfer files to remote computer plugin by Lee Katz";

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies tempdir=s help inbox=s debug force force-transfer numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      rsync     => 'rsync --version | grep version',
      ssh       => 'ssh -V 2>&1',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  # Remote pid file, to block multiple transfers
  $$settings{transfer_destination_string} =~ m/(.+?\@)?(.+?)\:(.+)?/;
  my($username,$url,$remotePath)=($1,$2,$3);
  if(!$url){
    logmsg "WARNING: I could not parse $$settings{transfer_destination_string} for a remote server!";
  }
  $username//="";
  $username=~s/\@$//;

  my $remotePid="$remotePath/.SneakerNet/pid.txt";

  # Check for the remote pid file
  my $pid=`ssh -q $username\@$url cat $remotePid 2>/dev/null`;
  logmsg "WARNING: I could not check for the remote pid file (pid: $pid): $!" if $?;
  $pid||=0;
  $pid+=0;
  
  my $numTries=0;
  while($pid > 0 && !$$settings{force}){
    logmsg "ERROR: there is either already a transfer in progress into target folder $remotePath or a previous iteration died.  The local pid is/was $pid. Run this script with --force to ignore this error.";
    logmsg "I will sleep 1 minute to try again.  Delete $remotePid on remote computer to avoid this warning.";
    sleep 60;
    $pid=`ssh -q $username\@$url cat $remotePid 2>/dev/null`;

    $numTries++;
    if($numTries > 10){
      logmsg "Gave up after $numTries tries. I'll continue the transfer anyway, and I'll remove $url:$remotePid.";
      command("ssh -q $username\@$url rm -vf $remotePid");
    }
  }

  # Make the pid file
  command("ssh -q $username\@$url 'mkdir -pv $remotePath/.SneakerNet; echo $$ > $remotePid'");

  transferFilesToRemoteComputers($dir,$settings);

  # Remove the remote pid file
  command("ssh -q $username\@$url rm -vf $remotePid");

  recordProperties($dir,{
    version=>$VERSION,
    dateSent=>strftime("%Y-%m-%d", localtime()),
    timeSent=>strftime("%H:%M:%S", localtime()),
  });

  return 0;
}

sub transferFilesToRemoteComputers{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  # Which files should be skipped according to Q/C?
  my $passfail=passfail($dir,$settings);
  
  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  SAMPLE:
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    # Find the taxon
    my $taxon=$$s{species} || $$s{taxon} || 'NOT LISTED';
    logmsg "The taxon of $sampleName is $taxon";
    my @route = (ref($$s{route}) eq 'ARRAY')?@{$$s{route}}:($$s{route});

    # If we are transferring...
    # if($$settings{'force-transfer'} || (scalar(@route) > 1 && grep {/calcengine/i} @route) ){
    # addressing issue #58
    if($$settings{'force-transfer'} || grep {/calcengine/i} @route ){
      # if this sample fails at all, then NEXT!
      for my $reason(keys(%{ $$passfail{$sampleName} })){
        if($$passfail{$sampleName}{$reason} == 1){
          logmsg "Failed $sampleName because failed category $reason";
          next SAMPLE;
        }
      }

      # Where does this sample get transferred?
      my $subfolder=$$s{taxonRules}{dest_subfolder};
      if(!defined $$s{taxonRules}{dest_subfolder}){
        if($$s{catchall_subfolder}){
          $subfolder=$$s{catchall_subfolder};
        } else {
          $subfolder="SneakerNet";
        }
      }

      # Add on these fastq files for transfer
      for my $fastq(@{ $$s{fastq} }){
        $filesToTransfer{$subfolder} .= "$fastq ";
      }

      logmsg "One route for sample $sampleName is the Calculation Engine";
    } else {
      logmsg "Note: The route for $sampleName was not listed in the sample sheet.";
    }
  }

  #die "ERROR: no files to transfer" if (!$filesToTransfer);
  logmsg "WARNING: no files will be transferred. Use --force-transfer to override." if(!keys(%filesToTransfer));

  # Make the transfers based on taxon.
  while(my($subfolder,$fileString)=each(%filesToTransfer)){

    logmsg "Transferring to $subfolder:\n  $fileString";
    next if($$settings{debug});
    eval{
      #print "rsync -q --no-motd --update -av --no-g $fileString $$settings{transfer_destination_string}/$subfolder/\n";
      command("rsync -av -q --no-motd --update -av --no-g --copy-links --chmod=Du=rwx,Dg=rx,Do=rx,Fu=rw,Fg=r,Fo=r $fileString $$settings{transfer_destination_string}/$subfolder/");
    };
    if($@){
      logmsg "ERROR: I could not transfer these files. If it is a permissions issue, one cause is if the destination file already exists but under a different username.";
    }
  }
}

sub usage{
  print "Find all reads directories under the inbox
  Usage: $0 MiSeq_run_dir
  --debug            No files will actually be transferred
  --force            Ignore some warnings
  --force-transfer   Transfer the reads despite the routing
                     entry in the spreadsheet.
  --version
";
  exit(0);
}

