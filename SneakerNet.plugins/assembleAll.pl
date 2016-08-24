#!/usr/bin/env perl
# Transfers files to a remote destination and QCs them beforehand.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use FindBin;
use lib "$FindBin::RealBin/../lib";
use SneakerNet qw/readConfig samplesheetInfo command logmsg fullPathToExec/;

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help tempdir=s debug numcpus=i)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  # Check for required executables
  for (qw(megahit prodigal run_assembly_filterContigs.pl run_prediction_metrics.pl)){
    fullPathToExec($_);
  }
 
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/assemblies";
  my $metricsOut=assembleAll($dir,$settings);

  return 0;
}

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo("$dir/SampleSheet.csv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");

    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly="$outdir/$sample.megahit.fasta";
    my $outgbk="$outdir/$sample.megahit.gbk";

    # Run the assembly
    if(!-e $outassembly){
      my $assembly=assembleSample($sample,$info,$settings);
      next if(!$assembly);

      # Save the assembly
      mkdir $outdir;
      command("run_assembly_filterContigs.pl -l 500 $assembly > $outassembly");
    }

    # Genome annotation
    if(!-e $outgbk){
      my $gbk=annotateFasta($sample,$outassembly,$settings);
      cp($gbk,$outgbk) or die "ERROR: could not copy $gbk to $outgbk: $!";
    }
  }
  
  # run assembly metrics with min contig size=0.5kb
  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  logmsg "Running metrics on the genbank files at $metricsOut";
  command("run_prediction_metrics.pl $dir/SneakerNet/assemblies/*/*.megahit.gbk > $metricsOut");
  
  return $metricsOut;
}

sub assembleSample{
  my($sample,$sampleInfo,$settings)=@_;

  my $R1=$$sampleInfo{fastq}[0];
  my $R2=$$sampleInfo{fastq}[1];
  if(!$R1){
    logmsg "Could not find R1 for $sample. Skipping.";
    return "";
  }
  if(!$R2){
    logmsg "Could not find R2 for $sample. Skipping.";
    return "";
  }

  logmsg "Assembling $sample";

  my $outdir="$$settings{tempdir}/$sample";
  
  system("rm -rf '$outdir'"); # make sure any previous runs are gone
  my $numcpus=$$settings{numcpus};
  $numcpus=2 if($numcpus < 2); # megahit requires at least two
  command("megahit -1 $R1 -2 $R2 --out-dir '$outdir' -t $numcpus 1>&2");
  die "ERROR with running megahit on $sample: $!" if $?;

  return "$outdir/final.contigs.fa";
}

# I _would_ use prokka, except it depends on having an up to date tbl2asn
# which is not really necessary for what I'm doing here.
sub annotateFasta{
  my($sample,$assembly,$settings)=@_;

  my $outdir="$$settings{tempdir}/$sample/prokka";
  system("rm -rf $outdir");
  mkdir $outdir;
  my $outgff="$outdir/prodigal.gff";
  my $outgbk="$outdir/prodigal.gbk";

  logmsg "Predicting genes on $sample with Prodigal";
  command("prodigal -q -i $assembly -o $outgff -f gff -g 11 1>&2");

  # Read the assembly sequence
  my %seqObj;
  my $seqin=Bio::SeqIO->new(-file=>$assembly);
  while(my $seq=$seqin->next_seq){
    $seqObj{$seq->id}=$seq;
  }
  $seqin->close;

  # Add seq features
  my $gffin=Bio::FeatureIO->new(-file=>$outgff);
  while(my $feat=$gffin->next_feature){
    # put the features onto the seqobj and write it to file
    my $id=$feat->seq_id;
    $seqObj{$id}->add_SeqFeature($feat);
  }
  $gffin->close;

  # Convert to gbk
  my $gbkObj=Bio::SeqIO->new(-file=>">$outgbk",-format=>"genbank");
  for my $seq(values(%seqObj)){
    $gbkObj->write_seq($seq);
  }
  $gbkObj->close;

  return $outgbk;
}

sub usage{
  "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  "
}

