#!/usr/bin/env perl
# Runs GP60 on a genome assembly

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir tempfile/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.0";
our $CITATION= "Gp60 plugin by Lee Katz. Uses gp60 by Alyssa Kelley.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      'blastn (BLAST+)'    => 'blastn -version | head -n 1',
      rm        => 'rm --version | head -n 1',
      'countGP60repeats.pl (GP60_Counter - private repo)' => 'echo GP60 counter unknown version',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/gp60";
  mkdir "$dir/SneakerNet/forEmail";

  my $db = $$settings{gp60db} || die "ERROR: need to add gp60db to settings.conf";

  my $reportHash = gp60ForAll($dir,$db,$settings);
  my $reportFile = "$dir/SneakerNet/forEmail/gp60.tsv";

  # Make a report for email
  open(my $outFh, '>', $reportFile) or die "ERROR: could not write to $reportFile: $!";
  print $outFh join("\t", qw(sample gp60))."\n";
  while(my($sample, $report) = each(%$reportHash)){
    open(my $inFh, '<', $report) or die "ERROR: could not read from $report: $!";
    my $result = <$inFh>;
    close $inFh;
    chomp($result);
    my @F = split(/\t/, $result);

    print $outFh join("\t", $sample, $F[3])."\n";
  }
  close $outFh;

  logmsg "Results found in $reportFile";

  recordProperties($dir,{version=>$VERSION, table=>$reportFile});

  return 0;
}

sub gp60ForAll{
  my($dir, $db, $settings) = @_;

  my @queueBuffer = ();

  # Find information about each genome
  logmsg "Reading sample tsv at $dir/samples.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    
    my $assembly=(glob("$dir/SneakerNet/assemblies/$sample/$sample.*.fasta"))[0];
    if(!$assembly){
      logmsg "ERROR: no assembly was found for $sample. Will not run analysis on this one sample.";
      next;
    }

    push(@queueBuffer, {assembly=>$assembly, sample=>$sample});
  }

  # make sure that the gp60 database has been formatted
  for my $dbFile("$db.nhr", "$db.nin", "$db.nsq"){
    if(! -e $dbFile){
      logmsg "WARNING: database was given as $db but it might not be formatted. Could not find file $dbFile";
      logmsg "WARNING: formatting the database myself.";
      command("makeblastdb -in $db -dbtype nucl");
      last;
      #die "ERROR: database was given as $db but it might not be formatted. Could not find file $dbFile";
    }
  }

  # Kick off the threads with the array buffer
  my $Q=Thread::Queue->new(@queueBuffer);
  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i]=threads->new(\&gp60Worker,$Q,$dir,$db,$settings);
    # also send a terminating signal
    $Q->enqueue(undef);
  }

  my %reportHash =();
  for(@thr){
    my $reports = $_->join;
    %reportHash = (%reportHash, %$reports);
  }

  return \%reportHash;
}

sub gp60Worker{
  my($Q, $dir, $db, $settings) = @_;

  my %report;

  while(defined(my $info = $Q->dequeue)){
    my $asm    = $$info{assembly};
    my $sample = $$info{sample};
    my($blastoutFh, $blastout) = tempfile("XXXXXX", SUFFIX=>".bls.tsv",  DIR=>$$settings{tempdir}, CLEANUP=>1);
    my($gp60Fh, $gp60)         = tempfile("XXXXXX", SUFFIX=>".gp60.tsv", DIR=>$$settings{tempdir}, CLEANUP=>1);
    
    command("blastn -query $asm -db $db -outfmt '6 qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send qlen slen qseq' | head -n 1 > $blastout");

    command("countGP60repeats.pl $blastout > $gp60");
    
    $report{$sample} = $gp60;
  }

  return \%report;
}

sub usage{
  print "Run gp60 analysis on a set of assemblies
  Usage: $0 MiSeq_run_dir
  --numcpus 1
";
  exit(0);
}

