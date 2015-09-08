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

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug)) or die $!;
  die usage() if($$settings{help} || !@ARGV);

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
    my $taxon=(keys(%{ $$s{species} }))[0];
    if($$s{route}{calcengine}){
      $filesToTransfer{$taxon}.=join(" ",glob("$dir/$sampleName*.fastq.gz"))." ";
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
    command("rsync --update -av $fileString edlb-sneakernet\@biolinux.biotech.cdc.gov:/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/out/Calculation_Engine/$subfolder/");
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

  # Try to associate samples to files
  while(my($samplename,$sampleinfo)=each(%sample)){
    my @possibleFastq=glob(dirname($samplesheet)."/$samplename*.fastq.gz");
    $sample{$samplename}{fastq}=\@possibleFastq;
  }

  return \%sample;
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
