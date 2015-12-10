#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;

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

  # Which files should be transferred?
  my %filesToTransfer=(); # hash keys are species names
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases
    my $taxon=$$s{species} || 'NOT LISTED';
    logmsg "The taxon of $sampleName is $taxon";
    if(grep {/calcengine/i} @{ $$s{route} }){
      $filesToTransfer{$taxon}.=join(" ",@{ $$s{fastq} })." ";
      logmsg "One route for sample $sampleName is the Calculation Engine";
    } else {
      logmsg "Note: The route for $sampleName was not listed in the sample sheet.";
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
    if($taxon =~ /Listeria|^L\.$/i){
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
    command("rsync --update -av $fileString $$settings{transfer_destination_string}/$subfolder/");
  }
}


sub usage{
  "Find all reads directories under the inbox
  Usage: $0 MiSeq_run_dir
  "
}

