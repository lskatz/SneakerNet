#!/usr/bin/env perl
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/move copy/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

use threads;
use Thread::Queue;

use lib "$FindBin::RealBin/../lib/perl5";
use List::MoreUtils qw/uniq/;
use SneakerNet qw/readMetrics exitOnSomeSneakernetOptions recordProperties readConfig logmsg samplesheetInfo_tsv command/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";
our $VERSION = "2.0";
our $CITATION = "Add read metrics by Lee Katz. Uses read metrics script in CG-Pipeline.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help force citation check-dependencies inbox=s debug test numcpus=i tempdir=s version)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir($0.".XXXXXX", TMPDIR=>1, CLEANUP=>1);
  logmsg "Tempdir is $$settings{tempdir}";
  
  my $dir=$ARGV[0];

  addReadMetrics($dir,$settings);

  # Mark this file as something to attach for an email later
  mkdir "$dir/SneakerNet" if(!-d "$dir/SneakerNet");
  mkdir "$dir/SneakerNet/forEmail" if(!-d "$dir/SneakerNet/forEmail");
  if(! -f "$dir/SneakerNet/forEmail/readMetrics.tsv"){
    link("$dir/readMetrics.tsv","$dir/SneakerNet/forEmail/readMetrics.tsv")
      or die "ERROR: could not hard link readMetrics.tsv to $dir/SneakerNet/forEmail/ - $!";
  }

  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{version=>$VERSION, table=>"$dir/SneakerNet/forEmail/readMetrics.tsv", mqc=>$rawMultiQC});

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/readMetrics.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/readMetrics_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Read metrics\"\n";
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

sub addReadMetrics{
  my($dir,$settings)=@_;

  return if(-e "$dir/readMetrics.tsv" && !$$settings{force});

  logmsg "Reading sample $dir/SampleSheet.csv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  logmsg "Running fast read metrics";

  my $Q=Thread::Queue->new(glob("$dir/*.fastq.gz"));
  my @thr;
  for (0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&readMetricsWorker, $Q, $settings);
    $Q->enqueue(undef);
  }

  for(@thr){
    $_->join;
  }

  system("cat $$settings{tempdir}/*/readMetrics.tsv | head -n 1 > $dir/readMetrics.tsv.tmp"); # header
  system("sort -k3,3n $$settings{tempdir}/*/readMetrics.tsv | uniq -u >> $dir/readMetrics.tsv.tmp"); # content
  die if $?;

  # edit read metrics to include genome sizes
  logmsg "Backfilling values in $dir/readMetrics.tsv";
  my $newReadMetrics;
  open(READMETRICS,"$dir/readMetrics.tsv.tmp") or die "ERROR: could not open $dir/readMetrics.tsv.tmp because $!";
  open(READMETRICSFINAL,">","$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv for writing: $!";

  # get the header and also put it into the final output file
  my $header=<READMETRICS>;
  print READMETRICSFINAL $header;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<READMETRICS>){
    chomp;
    # read in each line into the appropriate header
    my %h;
    @h{@header}=split(/\t/);

    # find the genome size based on the filename
    my $coverage=calculateCoverage(\%h,$sampleInfo,$settings);
    $h{coverage}=$coverage;
    $h{File}=basename($h{File});

    for(@header){
      print READMETRICSFINAL "$h{$_}\t";
    }
    print READMETRICSFINAL "\n";
  }
  close READMETRICSFINAL;
  close READMETRICS;

  # Clean up by removing the temporary file
  unlink("$dir/readMetrics.tsv.tmp");
}

sub readMetricsWorker{
  my($Q, $settings)=@_;
  
  my %metrics;

  my $tempdir=tempdir("worker.XXXXXX", DIR=>$$settings{tempdir}, CLEANUP=>1);
  while(defined(my $fastq=$Q->dequeue)){
    logmsg "read metrics for $fastq";
    my $metricsHashRef = readMetrics([$fastq]);
    %metrics = (%metrics, %$metricsHashRef);
  }

  # Write to the output file
  open(my $fh, ">>", "$tempdir/readMetrics.tsv") or die "ERROR: could not append to $tempdir/readMetrics.tsv: $!";
  print $fh join("\t", qw(File avgReadLength totalBases minReadLength maxReadLength avgQuality numReads coverage))."\n";
  while(my($fastq,$m) = each(%metrics)){
    print $fh join("\t", basename($fastq), $$m{avgReadLength}, $$m{totalBases}, $$m{minReadLength}, 
                         $$m{maxReadLength}, $$m{avgQuality}, $$m{numReads}, $$m{coverage}) . "\n";
  }
  close $fh;
}


# Use a hash of a line from a readmetrics output 
# to determine genome coverage.
sub calculateCoverage{
  my($h,$sampleInfo,$settings)=@_;

  # $h contains read metrics for this one row.

  my $file=basename($$h{File});
  my $samplename=$file || "";

  # If we don't have a sample name then it's probably because
  # the "sample name" is a filename.
  if(!$$sampleInfo{$samplename}){
    SAMPLENAME:
    for my $sn(keys(%$sampleInfo)){
      for my $fastq(@{ $$sampleInfo{$sn}{fastq} }){
        #logmsg basename($fastq)." eq $file";
        if(basename($fastq) eq $file){
          #logmsg "FOUND: $sn";
          $samplename = $sn;
          last SAMPLENAME;
        }
      }
    }
  }
  # If we still don't have this sample, then quit
  if(!$$sampleInfo{$samplename}){
    return 0;
  }

  # Find out if this file has an expected genome size from the Sample Sheet.
  my $expectedGenomeSize = $$sampleInfo{$samplename}{expectedgenomesize} || $$sampleInfo{$samplename}{taxonRules}{genomesize} || 0;
  if($expectedGenomeSize < 100){
    $expectedGenomeSize *= 10**6;
  }
  my $organism = $$sampleInfo{$samplename}{taxon};
  #die Dumper $sampleInfo, $organism, $samplename, $file if($samplename =~ /43410/);

  my $coverage=$$h{coverage} || 0; 

  # Recalculate coverage, if it's possible
  if($expectedGenomeSize > 0){
    $coverage=$$h{totalBases}/$expectedGenomeSize;
    $coverage=sprintf("%0.2f",$coverage); # round it
    logmsg "Decided that $$h{File} is $organism with expected genome size $expectedGenomeSize. Calculated coverage: $coverage";
  } else {
    logmsg "Warning: could not understand what organism $$h{File} belongs to. I tried to look it up by $samplename. Coverage was not recalculated.";
    #logmsg Dumper $sampleInfo, $file;
  }
  #die Dumper $$sampleInfo{$samplename}, $expectedGenomeSize, $$h{coverage}, $coverage;
  return $coverage;
}
 
################
# Utility subs #
################

sub usage{
  print "Find all reads directories under the inbox
  Usage: $0 runDir
  --debug # Show debugging information
  --numcpus  1
  --tempdir ''
  --force
  --version
  ";

  exit(0);
}

