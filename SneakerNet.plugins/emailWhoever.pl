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

use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig command logmsg/;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help numcpus=i inbox=s debug)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  emailWhoever($dir,$settings);

  return 0;
}

sub emailWhoever{
  my($dir,$settings)=@_;

  my $subdir=basename($dir);
  my $machineName=basename(dirname($dir));
  my $readMetrics="$dir/readMetrics.tsv";

  # Figure out who we are mailing
  # Start off with what emails are in the config
  my @to=flatten(@{ $$settings{emails} });

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
    my $body ="Please open the following attachments for QC information on $subdir.\n";
       $body.=" - TSV files can be opened in Excel\n";
       $body.=" - LOG files can be opened in Wordpad\n";
       $body.=" - HTML files can be opened in Internet Explorer\n";
       $body.="\nThis message was brought to you by SneakerNet!\n";

    my $email=Email::Stuffer->from($from)
                               ->subject($subject)
                               ->to($to)
                               ->text_body($body);

    for my $file(glob("$dir/SneakerNet/forEmail/*")){
      $email->attach_file($file);
    }

    my $was_sent=$email->send;

    if(!$was_sent){
      logmsg "Warning: Email was not sent to $to!";
    }
  }

  return \@to;
}


################
# Utility subs #
################

# http://stackoverflow.com/a/20359734
sub flatten {
  map { ref $_ ? flatten(@{$_}) : $_ } @_;
}

sub usage{
  "Find all reads directories under the inbox
  Usage: $0 [-i inboxDir/]
  -i dir  # choose a different 'inbox' to look at
  --debug # Show debugging information
  "
}
