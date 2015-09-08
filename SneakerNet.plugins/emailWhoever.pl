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
use Email::Stuffer;
use List::MoreUtils qw/uniq/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help inbox=s debug)) or die $!;
  die usage() if($$settings{help} || !@ARGV);

  my $dir=$ARGV[0];

  emailWhoever($dir,$settings);

  return 0;
}

sub emailWhoever{
  my($dir,$settings)=@_;

  my $subdir=basename($dir);
  my $readMetrics="$dir/readMetrics.tsv";

  # Figure out who we are mailing
  # Send here by default
  my @to=("gzu2\@cdc.gov","wwm8\@cdc.gov","pfge\@cdc.gov","wvt2\@cdc.gov","fid4\@cdc.gov");
  # Read the sample sheet for something like
  #                                   Investigator Name,ALS (IAU3)
  # And then make IAU3 into a CDC email.
  my $pocLine=`grep -m 1 'Investigator' $dir/SampleSheet.csv`;
  if($pocLine=~/\((.+)\)/){
    my $cdcids=$1;
    $cdcids=~s/\s+//g;
    my @cdcids=map {"$_\@cdc.gov"} split(/,/,$cdcids);
    push(@to,@cdcids);
  } else {
    logmsg "WARNING: could not parse the investigator line so that I could find CDC IDs";
  }

  # Send one email per recipient.
  for my $to(uniq(@to)){
    logmsg "To: $to";
    my $from="sequencermaster\@monolith0.edlb.cdc.gov";
    my $subject="$subdir QC";
    my $body ="Please open the following attachment in Excel for read metrics for run $subdir.\n";
       $body.="\nThis message was brought to you by SneakerNet!";

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
  --debug # Show debugging information
  "
}
