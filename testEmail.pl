use Email::Stuffer;
use strict;
use warnings;

sub logmsg {print "@_\n"}

my $settings={};
my $info={dir=>"/mnt/monolith0Data/RawSequenceData/M347/M347-15-020",subdir=>"M347-15-020"};
emailWhoever($info,$settings);

sub emailWhoever{
  my($info,$settings)=@_;

  my $subdir=$$info{subdir};
  my $readMetrics=$$info{dir}."/readMetrics.txt";

  # Figure out who we are mailing
  my @to=("gzu2\@cdc.gov");
  my $pocLine=`grep -m 1 'Investigator' $$info{dir}/SampleSheet.csv`;
  if($pocLine=~/\((.+)\)/){
    my $cdcids=$1;
    $cdcids=~s/\s+//g;
    my @cdcids=map {"$_\@cdc.gov"} split(/,/,$cdcids);
    push(@to,@cdcids);
  } else {
    logmsg "WARNING: could not parse the investigator line so that I could find CDC IDs";
  }

  # Send one email per recipient.
  for my $to(@to){
    logmsg "To: $to";
    my $from="sequencermaster\@monolith0.edlb.cdc.gov";
    my $subject="$subdir QC";
    my $body="Please open the following attachment in Excel for read metrics for run $subdir.\n";

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

