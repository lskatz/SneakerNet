#!/usr/bin/env perl
# Performs 7-gene MLST on all samples

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use threads;
use Thread::Queue;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "2.0";
our $CITATION= "MLST plugin by Lee Katz. Uses mlst by Torsten Seemann.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      mlst      => 'mlst --version',
      'blastn (BLAST+)'    => 'blastn -version | head -n 1',
      rm        => 'rm --version | head -n 1',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/mlst";
  mkdir "$dir/SneakerNet/forEmail";
  my $reportArr=mlst($dir,$settings);

  # Make a report for email
  my $mlstAttachment="$dir/SneakerNet/forEmail/mlst.tsv";
  open(my $fh, ">", $mlstAttachment) or die "ERROR: could not write to $mlstAttachment: $!";
  print $fh join("\t", qw(File mlst_scheme 7-gene_ST locus1 locus2 locus3 locus4 locus5 locus6 locus7))."\n";
  # Read each report and write to the email attachment
  for my $report(@$reportArr){
    open(my $reportFh, $report) or die "ERROR: could not read $report: $!";
    while(my $reportLine=<$reportFh>){
      print $fh $reportLine;
    }
    close $reportFh;
  }

  # Give some definitions for mlst notations
  print $fh "# n indicates an exact allele\n";
  print $fh "# ~n indicates a novel allele similar to n\n";
  print $fh "# n? indicates a partial match to known allele n\n";
  print $fh "# n,m indicates multiple alleles\n";
  print $fh "# - indicates that an allele is missing\n";
  print $fh "# For more details: https://github.com/tseemann/mlst\n";
  close $fh;

  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{version=>$VERSION, table=>"$dir/SneakerNet/forEmail/mlst.tsv", mqc=>$rawMultiQC});

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/mlst.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/mlst_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"MLST (7-gene)\"\n";
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

sub mlst{
  my($dir,$settings)=@_;

  my @mlstQueueBuffer=();
  
  # Find information about each genome
  logmsg "Reading sample tsv at $dir/samples.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    
    my $assembly=(glob("$dir/SneakerNet/assemblies/$sample/$sample.*.fasta"))[0];
    if(!$assembly){
      logmsg "ERROR: no assembly was found for $sample";
      next;
    }
    #my $assembly="$dir/SneakerNet/assemblies/$sample/$sample.megahit.fasta";
    push(@mlstQueueBuffer, {assembly=>$assembly, sample=>$sample});
  }

  # Kick off the threads with the array buffer
  my $Q=Thread::Queue->new(@mlstQueueBuffer);
  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i]=threads->new(\&mlstWorker,$Q,$dir,$settings);
    # also send a terminating signal
    $Q->enqueue(undef);
  }

  my @reportArr=();
  for(@thr){
    my $reportList=$_->join;
    push(@reportArr, @$reportList);
  }

  return \@reportArr;
}

sub mlstWorker{
  my($Q,$dir,$settings)=@_;
  my @reportArr=();

  while(defined(my $shortInfo=$Q->dequeue)){
    my $outdir="$dir/SneakerNet/mlst/$$shortInfo{sample}";
    my $outfile="$outdir/mlst.tsv";
    push(@reportArr,$outfile);

    # If we are forcing it, then remove the whole
    # output directory if it exists.
    if($$settings{force}){
      system("rm -rvf $outdir");
    }
    # Make the output directory
    mkdir $outdir;
    
    logmsg "MLST for $$shortInfo{sample}";
    if(!-e $outfile){
      my $mlst=mlstSample($$shortInfo{assembly},$settings);
      
      # Explicitly print the results instead of using the
      # shell syntax ">" because we are running `mlst`
      # based on whether or not the output file exists.
      open(my $outFh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
      print $outFh $mlst;
      close $outFh;
    }
  }
  return \@reportArr;
}

sub mlstSample{
  my($assembly,$settings)=@_;
  
  my $mlst=`mlst $assembly --nopath`;
  die if $?;
  return $mlst;
}
sub usage{
  print "Run MLST on all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --force        To overwrite previous results
  --version
";
  exit(0);
}

