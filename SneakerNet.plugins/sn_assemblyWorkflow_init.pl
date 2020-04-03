#!/usr/bin/env perl
# Initiates the assembly-only workflow

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

our $VERSION = "1.0";
our $CITATION = "Assembly workflow init by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      'prodigal'                      => "prodigal -v 2>&1 | grep -i '^Prodigal V'",
      'run_assembly_metrics.pl (CG-Pipeline)'       => "echo CG Pipeline version unknown",
      'run_prediction_metrics.pl (CG-Pipeline)'     => "echo CG Pipeline version unknown",
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];

  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/assemblies";
  mkdir "$dir/SneakerNet/forEmail";

  my $metricsOut = "$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  evaluateAssemblies($dir,$settings);
  recordProperties($dir,{version=>$VERSION, table=>$metricsOut});

  logmsg "Metrics can be found in $metricsOut";

  return 0;
}

sub evaluateAssemblies{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly="$outdir/$sample.init.fasta";
    my $outgbk="$outdir/$sample.init.gbk";

    # TODO if a gbk is given, flatten it out to an assembly,
    # and then for consistency, run prodigal on it.

    # Cp the assembly
    if(!-e $outassembly){
      mkdir $outdir;

      $$info{fastq}[0] or die "ERROR: no assembly was defined for $sample";
      my $inassembly = File::Spec->rel2abs($$info{fastq}[0], $dir);
      logmsg "Copying $inassembly => $outassembly";
      cp($inassembly, $outassembly) or die "ERROR copying $inassembly => $outassembly : $!";
    }

    # Genome annotation
    if(!-e $outgbk){
      mkdir "$outdir/prodigal"; # just make this directory right away
      my $gbk=annotateFasta($sample,$outassembly,$settings);
      cp($gbk,$outgbk) or die "ERROR: could not copy $gbk to $outgbk: $!";
    }
  }
  
  # run assembly metrics with min contig size=0.5kb
  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  if(!-e $metricsOut){

    logmsg "Running metrics on the genbank files at $metricsOut";

    my @thr;
    my $Q=Thread::Queue->new(glob("$dir/SneakerNet/assemblies/*/*.init.gbk"));
    for(0..$$settings{numcpus}-1){
      $thr[$_]=threads->new(\&predictionMetricsWorker,$Q,$settings);
      $Q->enqueue(undef);
    }
    for(@thr){
      $_->join;
    }
    
    command("echo $$settings{tempdir}/worker.*; ls $$settings{tempdir}/worker.*");
    command("cat $$settings{tempdir}/worker.*/metrics.tsv | head -n 1 > $metricsOut.tmp"); # header
    command("sort -k1,1 $$settings{tempdir}/worker.*/metrics.tsv | uniq -u >> $metricsOut.tmp"); # content
    mv("$metricsOut.tmp", $metricsOut);
  }
  
  return $metricsOut;
}

sub predictionMetricsWorker{
  my($Q,$settings)=@_;
  my $tempdir=tempdir("worker.XXXXXX", DIR=>$$settings{tempdir}, CLEANUP=>1);
  my $predictionOut="$tempdir/predictionMetrics.tsv";
  my $assemblyOut  ="$tempdir/assemblyMetrics.tsv";
  my $metricsOut   ="$tempdir/metrics.tsv";
  command("touch $predictionOut $assemblyOut");

  my $numMetrics=0;
  while(defined(my $gbk=$Q->dequeue)){
    # Metrics for the fasta: the fasta file has the same
    # base name but with a fasta extension
    my $fasta=$gbk;
    $fasta=~s/gbk$/fasta/;
    
    logmsg "gbk metrics for $gbk";
    command("run_prediction_metrics.pl $gbk >> $predictionOut");
    logmsg "asm metrics for $fasta";
    command("run_assembly_metrics.pl --allMetrics --numcpus 1 $fasta >> $assemblyOut");
    $numMetrics++;
  }
  # Don't do any combining if no metrics were performed
  return if($numMetrics==0);

  # Combine the files
  open(my $predFh,$predictionOut) or die "ERROR: could not read $predictionOut: $!";
  open(my $assemblyFh, $assemblyOut) or die "ERROR: could not read $assemblyOut: $!";
  open(my $metricsFh, ">", $metricsOut) or die "ERROR: could not write to $metricsOut: $!";

  # Get and paste the header into the output metrics file
  my $predHeader=<$predFh>;
  my $assemblyHeader=<$assemblyFh>;
  chomp($predHeader,$assemblyHeader);
  my @predHeader=split(/\t/,$predHeader);
  my @assemblyHeader=split(/\t/,$assemblyHeader);
  print $metricsFh "File\tgenomeLength";
  for(@predHeader, @assemblyHeader){
    next if($_ =~ /File|genomeLength/);
    print $metricsFh "\t".$_;
  }
  print $metricsFh "\n";

  # Get and paste the metrics from asm and pred into the output file
  while(my $predLine=<$predFh>){
    my $assemblyLine=<$assemblyFh>;
    chomp($predLine,$assemblyLine);
    
    my %F;
    @F{@predHeader}=split(/\t/,$predLine);
    @F{@assemblyHeader}=split(/\t/,$assemblyLine);
    # In case there was no assembly or predictions, add in a default value
    for(@predHeader,@assemblyHeader){
      $F{$_} //= 'NA';
    }
    # Also take care of genomeLength specifically
    $F{genomeLength} //= "NA";
    $F{File}=basename($F{File},qw(.gbk .fasta));

    print $metricsFh $F{File}."\t".$F{genomeLength};
    for(@predHeader, @assemblyHeader){
      next if($_ =~ /File|genomeLength/);
      print $metricsFh "\t".$F{$_}
    }
    print $metricsFh "\n";
  }
  close $_ for($predFh,$assemblyFh,$metricsFh);
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

  mkdir $outdir; # Skesa will need this container
  #command("skesa --cores $numcpus --gz --fastq $R1,$R2 > $$settings{tempdir}/$sample.fasta");
  # faster skesa command taken from 
  #  https://github.com/ncbi/SKESA/issues/11#issuecomment-429007711
  system("skesa --fastq $R1,$R2 --steps 1 --kmer 51 --cores $numcpus --vector_percent 1 > $$settings{tempdir}/$sample.fasta");
  if($?){
    command("skesa --fastq $R1,$R2 --cores $numcpus > $$settings{tempdir}/$sample.fasta");
  }
  return "$$settings{tempdir}/$sample.fasta";

  #$numcpus=2 if($numcpus < 2); # megahit requires at least two
  #command("megahit -1 $R1 -2 $R2 --out-dir '$outdir' -t $numcpus 1>&2");
  #die "ERROR with running megahit on $sample: $!" if $?;

  #return "$outdir/final.contigs.fa";
}

# I _would_ use prokka, except it depends on having an up to date tbl2asn
# which is not really necessary for what I'm doing here.
sub annotateFasta{
  my($sample,$assembly,$settings)=@_;

  my $outdir="$$settings{tempdir}/$sample/prodigal";
  system("rm -rf $outdir");
  mkdir "$$settings{tempdir}/$sample";
  mkdir $outdir;
  my $outgff="$outdir/prodigal.gff";
  my $outgbk="$outdir/prodigal.gbk";

  logmsg "Predicting genes on $sample with Prodigal";
  eval{
    command("prodigal -q -i $assembly -o $outgff -f gff -g 11 1>&2");
  };
  # If there is an issue, push ahead with a zero byte file
  if($@){
    logmsg "There was an issue with predicting genes on sample $sample. This might be caused by a poor assembly.";
    open(my $fh, ">", $outgff) or die "ERROR: could not write to $outgff: $!";
    close $fh;
  }

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
  print "Initiates the assembly-only workflow by placing the files in the right places and running gene prediction. Also creates assembly/annotation metrics.
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
  ";
  exit(0);
}

