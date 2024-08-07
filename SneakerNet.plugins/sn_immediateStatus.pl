#!/usr/bin/env perl
# Performs an immediate check on whether a run is good to go

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Spec;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

use Text::Levenshtein qw/distance/;

our $VERSION = "2.0";
our $CITATION= "Immediate status report by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version emails=s citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  my @exe = qw(sendmail);
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe       => \@exe,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/forEmail";

  my $errHash = doubleCheckRun($dir,$settings);

  # Record all warnings/errors into a file
  my $numErrors = 0;
  my $numWarnings = 0;
  my $outfile = "$dir/SneakerNet/forEmail/immediateReaction.tsv";
  open(my $fh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  print $fh join("\t", qw(ErrType Sample ErrKeyword Error))."\n";
  for my $errType(sort keys(%$errHash)){
    for my $sample(sort keys(%{ $$errHash{$errType} })){
      for my $errKeyword(sort keys(%{ $$errHash{$errType}{$sample} })){
        print $fh join("\t",
          $errType, $sample, $errKeyword, 
          $$errHash{$errType}{$sample}{$errKeyword}
        );
        print $fh "\n";

        if($errType eq 'fastq'){
          $numWarnings++;
        }
        elsif($errType eq 'sample'){
          $numErrors++;
        }

      }
    }
  }
  close $fh;

  # Record any errors that might have happened.
  # If no errors occured, use bool "0" and if some
  # errors, then make a message.
  my $errorMsg  = "0";
  my $warningMsg= "0";
  if($numErrors > 0){
    $errorMsg = "$numErrors errors";
  }
  if($numWarnings > 0){
    $warningMsg = "$numWarnings warnings";
  }

  my @to = ();
  if($$settings{emails}){
    my @tmp = split(/,/, $$settings{emails});
    push(@to, @tmp);
  }
  elsif(ref($$settings{'default.emails'}) eq 'ARRAY'){
    push(@to, @{ $$settings{'default.emails'} });
  } else {
    push(@to, $$settings{'default.emails'})
  }
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
    push(@to, @email);
  }
  my $to = join(",",@to);
  logmsg "Sending an email to @to";

  my $from=$$settings{from} || die "ERROR: need to set 'from' in the settings.conf file!";
  my $subject="Initial SneakerNet status for ".basename(File::Spec->rel2abs($dir));
  my $body = "If you see errors below, please contact the bioinformatics team with your run number and when you deposited the run. Include this file in your message.\nRun can be found at ".File::Spec->abs2rel($dir)."\n";
     $body.= "\nDocumentation can be found at https://github.com/lskatz/SneakerNet/blob/master/docs/plugins/sn_immediateStatus.pl.md\n";
  
  my $emailFile = "$$settings{tempdir}/email.txt";
  open(my $fh2, ">", $emailFile) or die "ERROR: could not write to $emailFile: $!";
  print $fh2 "To: $to\n";
  print $fh2 "From: $from\n";
  print $fh2 "Subject: $subject\n";
  print $fh2 "\n";
  print $fh2 "$body\n\n";
  append_attachment($fh2, $outfile);
  close $fh2;

  command("sendmail -t < $emailFile");

  recordProperties($dir,{
      version=>$VERSION, reportTo=>$to, 
      warnings=>$warningMsg, errors=>$errorMsg,
      table=>"$dir/SneakerNet/forEmail/immediateReaction.tsv",
      exe  => \@exe,
    });

  return 0;
}

sub doubleCheckRun{
  my($dir,$settings)=@_;

  my %errHash;

  if(!-e "$dir/samples.tsv"){
    $errHash{samplesSheetMissing} = 1;
  }

  my @fastq = sort glob("$dir/*.fastq.gz");
  my %fastqIndex;
  $fastqIndex{$_}=1 for(@fastq);

  logmsg "Reading sample tsv at $dir/samples.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  my @sample = sort keys(%$sampleInfo);
  for my $sample(@sample){
    my $info = $$sampleInfo{$sample};

    if(!$$info{fastq} || !ref($$info{fastq})){
      $errHash{sample}{$sample}{fastqNotFound}=1;
    } elsif(@{ $$info{fastq} } < 1){
      $errHash{sample}{$sample}{fastqNotFound}=1;
    } elsif(@{ $$info{fastq} } < 2){
      print Dumper $$info{fastq};
      $errHash{sample}{$sample}{pairMissing}="Sample $sample has only a single read: ".join(" ",@{ $$info{fastq} });
    }

    if($errHash{sample}{$sample}{fastqNotFound}){
      my $nearest = nearestString($sample, \@fastq, $settings);
      #my $tf = Text::Fuzzy->new($sample);
      #my $nearest = $tf->nearestv(\@fastq);
      $errHash{sample}{$sample}{fastqNotFound} = "For sample $sample, the closest named fastq files are $nearest";
    }

  }

  # now see if there are fastq files that are not mentioned
  for my $filename(@fastq){
    next if($filename=~/Undetermined/i);
    my $basename = basename($filename);
    $basename =~ s/_.*//;
    if(!$$sampleInfo{$basename}){
      #my $tf = Text::Fuzzy->new($filename);
      #my $nearest = $tf->nearestv(\@sample);
      my $nearest = nearestString($filename, \@sample, $settings);
      $errHash{fastq}{$filename}{sampleNotFound} = "Found fastq $filename but no entry in the sample sheet matching $basename.  Did you mean $nearest?";
    }
  }

  return \%errHash;
}

sub nearestString{
  my($query, $strings, $settings) = @_;

  my $minDist = ~0; # max int
  my $nearest = "";
  for my $string(@$strings){
    my $dist = distance($query, $string);
    if($dist < $minDist){
      $minDist = $dist;
      $nearest = $string;
    }
  }

  return $nearest;
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
  print "Double check a run and its completeness. Email a report.
  Usage: $0 MiSeq_run_dir
  --emails  ''   email1,[email2...]
  --version
";
  exit(0);
}

