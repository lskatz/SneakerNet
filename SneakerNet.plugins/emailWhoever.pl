#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Find qw/find/;
use Cwd qw/realpath/;
use File::Temp qw/tempdir/;
use MIME::Base64 qw/encode_base64/;
use POSIX qw/strftime/;
use IO::Compress::Zip qw(zip $ZipError);

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

use Config::Simple;
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig passfail command logmsg version/;
use List::MoreUtils qw/uniq/;

our $VERSION = "3.4";
our $CITATION= "Email whoever by Lee Katz";

my $snVersion=version();

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(citation check-dependencies version help force numcpus=i debug tempdir=s email-only|email|just=s)) or die $!;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  #$$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1);

  my @exe = qw(sendmail);
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe       => \@exe,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;

  my $dir=$ARGV[0];

  my $to = emailWhoever($dir,$settings);

  recordProperties($dir,{
    version=>$VERSION, 
    reportSentTo=>join(", ", @$to),
    dateSent=>strftime("%Y-%m-%d", localtime()),
    timeSent=>strftime("%H:%M:%S", localtime()),
    exe => \@exe,
  });

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/assemblyMetrics_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Assembly metrics\"\n";
  print $outFh "#description: \"$plugin v$VERSION $docLink $pluginLink\"\n";
  print $outFh "#anchor: '$anchor'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  while(<$fh>){
    next if(/^#/);
    print $outFh $_;
  }
  close $fh;

  return $outtable;
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
  if(-e "$dir/SampleSheet.csv"){
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
  } else {
    logmsg "WARNING: SampleSheet.csv not found; not parsing for additional emails. More than likely there are some email addresses anyway in snok.txt or emails.conf.";
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
  my $body ="Please see below for QC information on $runName.\n\n";
     $body.="For more details, please see the other attachments.\n";
     $body.=" - TSV files can be opened in Excel\n";
     $body.=" - LOG files can be opened in Wordpad, Notepad++, or VSCode\n";
     $body.=" - HTML files can be opened in Edge\n";
     $body.=" - Full path: ".realpath($dir)."/SneakerNet\n";
     $body.="\nThis message was brought to you by SneakerNet v$snVersion!\n";
     $body.="Documentation can be found at https://github.com/lskatz/SneakerNet\n";

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

  my $emailFile = "$$settings{tempdir}/email.txt";
  open(my $fh, ">", $emailFile) or die "ERROR: could not write to $emailFile: $!";
  print $fh "To: $to\n";
  print $fh "From: $from\n";
  print $fh "Subject: $subject\n";
  print $fh "\n";
  print $fh "$body\n";

  # Save a list of files to be attached
  my @attachment;
  push(@attachment,
    zip_directory("$dir/SneakerNet/forEmail", "$$settings{tempdir}/$runName.zip"),
    zip_file("$dir/SneakerNet/forEmail/report.html", $settings),
    zip_file("$dir/SneakerNet/forEmail/multiqc_report.html", $settings),
  );
  # Append log attachments to the email text file
  for my $file(glob("$dir/*.log")){
    next if(!-f $file);
    push(@attachment,
      zip_file($file, $settings),
    );
  }

  # Test attachments to make sure we want to keep them
  my @finalAttachment;
  for my $file(@attachment){
    if(-s $file > 1e7){
      logmsg "NOTE: $file is too big. I will not attach it.";
    } else {
      push(@finalAttachment, $file);
    }
  }

  if(scalar(@attachment) != scalar(@finalAttachment)){
    print $fh "\nSome files were not attached, probably because they are too large. \n";
    print $fh "You can view these files in " . realpath($dir) . " .\n";
    for my $file(@attachment){
      if(!grep{$_ eq $file} @finalAttachment){
        print $fh " - " . basename($file) . "\n";
      }
    }
  }

  # Finally, attach the files
  for my $file(@finalAttachment){
    append_attachment($fh, $file);
  }
  
  close $fh;
  command("sendmail -t < $emailFile");

  return \@to;
}


################
# Utility subs #
################

# http://stackoverflow.com/a/20359734
sub flatten {
  map { ref $_ ? flatten(@{$_}) : $_ } @_;
}

sub zip_directory {
    my($dir, $zip_file) = @_;

    my @files;

    # Find all files in the directory
    find(sub {
      # Skip if not a file
      return if(!-f);
      # Skip if it's an HTML (presumably a report file)
      return if($_ =~ /\.html$/);
      # Record the list of files in the forEmail directory
      push @files, $File::Find::name;
    }, $dir);

    # Zip the directory but exclude HTML files which are presumably report files
    #   -9 best compression
    #   -y store symlink and do not resolve
    #   -r recursive
    #   -q quiet
    system("cd $dir && zip --exclude '*.html' -9 -y -r -q $zip_file ./");
    die "Zip failed" if $?;

    return $zip_file;
}

sub zip_file {
  my($file, $settings) = @_;

  my $zip_file = $$settings{tempdir}."/".basename($file).".zip";

  # zip the file but exclude directory information
  my $dir = dirname($file);
  system("cd $dir && zip -9 -y -q -j $zip_file ". basename($file));
  die "Zip failed on $zip_file" if $?;

  #system("zip -9 -y -q $zip_file $file");
  #die "Zip failed on $zip_file" if $?;

  return $zip_file;
}

# Add an attachment to an email file handle
sub append_attachment {
    my ($fh, $file_path) = @_;

    # Encode the attachment content using base64 encoding
    my $attachment_name = basename($file_path);

    open(my $attachment_fh, "<", $file_path) or die "Failed to open attachment file $file_path: $!";
    binmode $attachment_fh;
    my $attachment_content = do { local $/; <$attachment_fh> };
    close $attachment_fh;

    my $encoded_content = pack("u", $attachment_content);
    die "Failed to encode attachment content from $file_path: $!" if $?;
    
    print $fh "begin 644 $attachment_name\n";
    print $fh $encoded_content . "\n";
    print $fh "end\n";

    # Print a newline to separate MIME parts
    print $fh "\n";
}


sub usage{
  print "Email a SneakerNet run's results
  Usage: $0 run-dir
  --debug          Show debugging information
  --numcpus     1  Number of CPUs (has no effect on this script)
  --email-only  '' Choose the email to send the report to instead
                   of what is supplied.
  --version
  ";
  exit(0);
}

