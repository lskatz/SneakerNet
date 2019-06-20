#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use Cwd qw/realpath/;
use File::Temp;
use FindBin;
use Config::Simple;
use POSIX qw/strftime/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readConfig passfail command logmsg version/;
use Email::Stuffer;
use List::MoreUtils qw/uniq/;

our $VERSION = "1.1";

my $snVersion=version();

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version help numcpus=i debug tempdir=s email-only=s)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  my $to = emailWhoever($dir,$settings);

  recordProperties($dir,{
    version=>$VERSION, 
    reportSentTo=>join(", ", @$to),
    dateSent=>strftime("%Y-%m-%d", localtime()),
    timeSent=>strftime("%H:%M:%S", localtime()),
  });

  return 0;
}

sub emailWhoever{
  my($dir,$settings)=@_;

  my $runName=basename(realpath($dir));
  my $machineName=basename(dirname($dir));
  my $readMetrics="$dir/readMetrics.tsv";

  # Figure out who we are mailing.
  # Start off with what emails are in the config, but pay
  # attention to whether there is just one or many.
  my @to;
  if(ref($$settings{'default.emails'}) eq "ARRAY"){
    push(@to, @{$$settings{'default.emails'}});
  } else {
    push(@to, $$settings{'default.emails'});
  }

  # Read the sample sheet for something like
  #                                   Investigator Name,ALS (IAU3)
  #     or:                           Investigator Name,ALS (IAU3;GZU2)
  # And then make IAU3 into a CDC email.
  my $pocLine=`grep -m 1 'Investigator' $dir/SampleSheet.csv`;
  if($pocLine=~/\((.+)\)/){
    my $cdcids=$1;
    $cdcids=~s/\s+//g; # remove whitespace
    for my $email(split(/[;,]/,$cdcids)){
      if($email !~ /\@/){
        $email.="\@cdc.gov";
      }
      push(@to,$email);
    }
  } else {
    logmsg "WARNING: could not parse the investigator line so that I could find CDC IDs";
  }

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
    push(@to,@email);
  }

  if($$settings{'email-only'}){
    @to=($$settings{'email-only'});
  }

  # Which files failed?
  my $failure=passfail($dir,$settings);

  # Get the email together for sending
  my $to=join(",",@to);
  logmsg "To: $to";
  my $from=$$settings{from} || die "ERROR: need to set 'from' in the settings.conf file!";
  my $subject="$runName QC";
  my $body ="Please see the report.html file for QC information on $runName.\n\n";
     $body.="For more details, please see the other attachments.\n";
     $body.=" - TSV files can be opened in Excel\n";
     $body.=" - LOG files can be opened in Wordpad\n";
     $body.=" - HTML files can be opened in Internet Explorer\n";
     $body.="\nThis message was brought to you by SneakerNet v$snVersion!\n";

  # Failure messages in the body
  $body.="\nAny samples that have failed QC as shown in passfail.tsv are listed below.\n";
  for my $fastq(keys(%$failure)){
    my $failureMessage="";
    for my $failureCategory(keys(%{$$failure{$fastq}})){
      if($$failure{$fastq}{$failureCategory} == 1){
        $failureMessage.=$fastq."\n";
        last; # just list a given failed fastq once
      }
    }
    $body.=$failureMessage;
  }

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
  "Email a SneakerNet run's results
  Usage: $0 run-dir
  --debug          Show debugging information
  --numcpus     1  Number of CPUs (has no effect on this script)
  --email-only  '' Choose the email to send the report to instead
                   of what is supplied.
  --version
  "
}
