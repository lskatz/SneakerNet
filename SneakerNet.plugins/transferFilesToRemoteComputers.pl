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

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig samplesheetInfo command logmsg passfail/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug force force-transfer numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
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
  
  if($pid > 0 && !$$settings{force}){
    die "ERROR: there is either already a transfer in progress into target folder $remotePath or a previous iteration died.  The local pid is/was $pid. Run this script with --force to ignore this error.";
  }

  # Make the pid file
  command("ssh -q $username\@$url 'mkdir $remotePath/.SneakerNet; echo $$ > $remotePid'");

  transferFilesToRemoteComputers($dir,$settings);

  # Remove the remote pid file
  command("ssh -q $username\@$url rm $remotePid");

  return 0;
}

sub transferFilesToRemoteComputers{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);

  # Which files should be skipped according to Q/C?
  my $passfail=passfail($dir,$settings);
  
  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases
    my $taxon=$$s{species} || 'NOT LISTED';
    logmsg "The taxon of $sampleName is $taxon";
    if($$settings{'force-transfer'} || grep {/calcengine/i} @{ $$s{route} }){
      FASTQ: for(@{ $$s{fastq} }){
        # Write out the status
        # The key of each filename is its basename
        my $f=basename($_);

        for my $reason(keys(%{ $$passfail{$f} })){
          if($$passfail{$f}{$reason} == 1){
            logmsg "Failed $f because $reason";
            next FASTQ;
          }
        }
        my $subfolder=$$s{taxonRules}{dest_subfolder} || "SneakerNet";
        $filesToTransfer{$subfolder}.=$_." ";
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
      command("rsync -q --no-motd --update -av --no-g $fileString $$settings{transfer_destination_string}/$subfolder/");
    };
    if($@){
      logmsg "ERROR: I could not transfer these files. If it is a permissions issue, one cause is if the destination file already exists but under a different username.";
    }
  }
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 MiSeq_run_dir
  --debug            No files will actually be transferred
  --force            Ignore some warnings
  --force-transfer   Transfer the reads despite the routing
                     entry in the spreadsheet.
  "
}

