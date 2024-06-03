#!/usr/bin/env perl
# Assemble genomes for Cryptosporidum

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

our $VERSION = "1.1";
our $CITATION= "Cryptosporidium assembly plugin by Lee Katz. Uses SHOvill.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  my @exe => (
    'run_prediction_metrics.pl', 'run_assembly_metrics.pl',
    qw(seqtk cat sort head uniq touch shovill prodigal),
    # shovill requires a ton of things:
    qw(seqtk pigz mash trimmomatic lighter flash
    spades.py skesa bwa samtools samclip java pilon),
    # not using megahit or velvet in this instance of shovill
    qw(megahit megahit_toolkit velveth velvetg),
  );
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe       => \@exe,
    }, $settings,
  );

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  logmsg "Temporary directory is at $$settings{tempdir}";

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";

  if($$settings{force}){
    logmsg "--force was given; removing any existing assembly results.";
    system("rm -rf $dir/SneakerNet/assemblies $dir/forEmail/assemblyMetrics.tsv");
  }
 
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/assemblies";
  my $metricsOut=assembleAll($dir,$settings);
  logmsg "Metrics can be found in $metricsOut";

  recordProperties($dir,{exe=>\@exe, version=>$VERSION,table=>$metricsOut});

  return 0;
}

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    logmsg "ASSEMBLE SAMPLE $sample";

    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly="$outdir/$sample.shovill.skesa.fasta";
    my $outgbk="$outdir/$sample.shovill.skesa.gbk";
    #my $outassembly="$outdir/$sample.megahit.fasta";
    #my $outgbk="$outdir/$sample.megahit.gbk";

    # Run the assembly
    if(!-e $outassembly){
      my $assembly=assembleSample($sample,$info,$settings);
      next if(!$assembly);
      next if(! -s $assembly);

      # Save the assembly
      mkdir $outdir;
      mkdir "$outdir/prodigal"; # just make this directory right away
      command("seqtk seq -L 500 $assembly > $outassembly");
    }

    # Genome annotation
    logmsg "PREDICT SAMPLE GENES $sample";
    if(!-e $outgbk){
      my $gbk=annotateFasta($sample,$outassembly,$settings);
      cp($gbk,$outgbk) or die "ERROR: could not copy $gbk to $outgbk: $!";
    }

  }
  
  # run assembly metrics with min contig size=0.5kb
  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  logmsg "Running metrics on the genbank files at $metricsOut";

  my @thr;
  my $Q=Thread::Queue->new(glob("$dir/SneakerNet/assemblies/*/*.shovill.skesa.gbk"));
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&predictionMetricsWorker,$Q,$settings);
    $Q->enqueue(undef);
  }
  for(@thr){
    $_->join;
  }
  
  command("cat $$settings{tempdir}/worker.*/metrics.tsv | head -n 1 > $metricsOut"); # header
  command("sort -k1,1 $$settings{tempdir}/worker.*/metrics.tsv | uniq -u >> $metricsOut"); # content
  
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

  my($R1,$R2) = @{ $$sampleInfo{fastq} };

  if(!$R1 || !-e $R1){
    logmsg "Could not find R1 for $sample. Skipping.";
    return "";
  }
  if(!$R2 || !-e $R2){
    logmsg "Could not find R2 for $sample. Skipping.";
    return "";
  }

  # Dealing with a small file size
  my $rawReadsFilesize = (stat($R1))[7] + (stat($R2))[7];
  if($rawReadsFilesize < 1000){
    if($$settings{force}){
      logmsg "Fastq files are tiny, but forcing because of --force.";
    } else {
      logmsg "Fastq files are too tiny. Skipping. Force with --force.";
      return "";
    }
  }

  logmsg "Assembling $sample";

  my $outdir="$$settings{tempdir}/$sample";
  
  system("rm -rf '$outdir'"); # make sure any previous runs are gone

  eval{
    command("shovill --outdir $outdir --R1 $R1 --R2 $R2 --assembler skesa --cpus $$settings{numcpus} skesa");
  };
  if($@){
    logmsg "shovill failed!\n$@";
    return "";
  }

  return "$outdir/contigs.fa";
}

# I _would_ use prokka, except it depends on having an up to date tbl2asn
# which is not really necessary for what I'm doing here.
sub annotateFasta{
  my($sample,$assembly,$settings)=@_;

  # Ensure a clean slate
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
  print "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
";
  exit(0);
}

