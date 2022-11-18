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

our $VERSION = "2.6.0";
our $CITATION= "Assembly plugin by Lee Katz. Uses SHOvill.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      'run_prediction_metrics.pl (CG-Pipeline)'     => "echo CG Pipeline version unknown",
      'run_assembly_metrics.pl (CG-Pipeline)'       => "echo CG Pipeline version unknown",
      cat                             => 'cat --version | head -n 1',
      sort                            => 'sort --version | head -n 1',
      head                            => 'head --version | head -n 1',
      uniq                            => 'uniq --version | head -n 1',
      touch                           => 'touch --version | head -n 1',
      shovill                         => 'shovill --version',
      prodigal                        => "prodigal -v 2>&1 | grep -i '^Prodigal V'",

      # shovill requires a ton of things:
      'seqtk'       => 'seqtk 2>&1 | grep Version',
      'pigz'        => 'pigz --version 2>&1',
      'mash'        => 'mash --version 2>&1',
      'pigz'        => 'pigz --version 2>&1',
      'trimmomatic' => 'trimmomatic -version 2>&1 | grep -v _JAVA',
      'lighter'     => 'lighter -v 2>&1',
      'flash'       => 'flash --version 2>&1 | grep FLASH',
      'spades.py'   => 'spades.py  --version 2>&1',
      'skesa'       => 'skesa --version 2>&1 | grep SKESA',
      'bwa'         => 'bwa 2>&1 | grep Version:',
      'samtools'    => 'samtools 2>&1 | grep Version:',
      'samclip'     => 'samclip --version 2>&1',
      'java'        => 'java -version 2>&1 | grep version',
      'pilon'       => 'pilon --version 2>&1 | grep -v _JAVA',

      # not using megahit or velvet in this instance of shovill
      'megahit'     => 'megahit --version 2>&1',
      'megahit_toolkit' => 'megahit_toolkit dumpversion 2>&1',
      'velveth'     => 'velveth 2>&1 | grep Version',
      'velvetg'     => 'velvetg 2>&1 | grep Version',
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

  recordProperties($dir,{version=>$VERSION,table=>$metricsOut,
    "Minimum contig length for assembly metrics" => "500bp",
  });

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
    my $outgraph="$outdir/$sample.shovill.skesa.gfa";
    #my $outassembly="$outdir/$sample.megahit.fasta";
    #my $outgbk="$outdir/$sample.megahit.gbk";

    # Run the assembly
    if(!-e $outassembly){
      my $assemblyDir = assembleSample($sample,$info,$settings);
      my $assembly    = "$assemblyDir/contigs.fa";
      next if(! -d $assemblyDir);
      next if(! -e $assembly);
      next if(! -s $assembly);

      # Save the assembly
      mkdir $outdir;
      mkdir "$outdir/shovill";
      for my $srcFile(glob("$assemblyDir/*")){
        my $target = "$outdir/shovill/".basename($srcFile);
        cp($srcFile, $target) or die "ERROR copying $srcFile => $target\n  $!";
      }
      cp($assembly, $outassembly) or die "ERROR copying $assembly => $outassembly\n  $!";

      mkdir "$outdir/prodigal"; # just make this directory right away
    }
    if(!-e $outgraph){
      my $graphdir = skesaToGraph($sample,$info,$outassembly,$settings);
      cp("$graphdir/contigs.gfa", $outgraph) or die "ERROR copying $graphdir/contigs.gfa => $outgraph\n  $!";
    }

    # Genome annotation
    logmsg "PREDICT GENES FOR SAMPLE $sample";
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
  #command("sort -k1,1 $$settings{tempdir}/worker.*/metrics.tsv | uniq -u >> $metricsOut"); # content
  command("tail -q -n +2 $$settings{tempdir}/worker.*/metrics.tsv | sort >> $metricsOut");
  
  return $metricsOut;
}

sub predictionMetricsWorker{
  my($Q,$settings)=@_;
  my $tempdir=tempdir("worker.XXXXXX", DIR=>$$settings{tempdir}, CLEANUP=>1);
  my $metricsOut    ="$tempdir/metrics.tsv";  # Combined metrics

  #command("touch $predictionOut $assemblyOut");

  my $numMetrics=0;
  while(defined(my $gbk=$Q->dequeue)){
    my $dir = dirname($gbk);
    my $predictionOut ="$dir/predictionMetrics.tsv";
    my $assemblyOut   ="$dir/assemblyMetrics.tsv";
    my $coverageOut   ="$dir/depth.tsv.gz";
    my $effectiveCoverage = "$dir/effectiveCoverage.tsv";

    # Metrics for the fasta: the fasta file has the same
    # base name but with a fasta extension
    my $fasta=$gbk;
    $fasta=~s/gbk$/fasta/;
    my $bam  =$gbk;
    $bam  =~s/gbk$/bam/;
    
    # Bring the shovill bam file over
    my $srcBam = dirname($fasta)."/shovill/shovill.bam";
    if(!-e $bam){
      logmsg "hard linking: $srcBam => $bam";
      link($srcBam, $bam) or die "ERROR: could not hard link $srcBam => $bam: $!";
    }
    
    logmsg "gbk metrics for $gbk";
    command("run_prediction_metrics.pl $gbk > $predictionOut");
    logmsg "asm metrics for $fasta";
    command("run_assembly_metrics.pl --allMetrics --numcpus 1 $fasta > $assemblyOut");
    logmsg "Coverage for $bam";
    command("samtools depth -aa $bam | gzip -c > $coverageOut");

    logmsg "Effective coverage for $bam";
    my $total = 0;
    my $numSites = 0;
    open(my $fh, "zcat $coverageOut |") or die "ERROR: could not zcat $coverageOut: $!";
    while(my $line = <$fh>){
      chomp($line);
      my(undef, undef, $depth) = split(/\t/, $line);
      $total+=$depth;
      $numSites++;
    }
    close $fh;
    open(my $outFh, ">", $effectiveCoverage) or die "ERROR: could not write to $effectiveCoverage: $!";
    print $outFh join("\t", qw(File effectiveCoverage))."\n";
    print $outFh join("\t", $bam, sprintf("%0.2f", $total/$numSites))."\n";
    close $outFh;

    # Combine the metrics into one file
    # We know that each file is one header and one values line
    my %metric;
    for my $file ($predictionOut, $effectiveCoverage, $assemblyOut){
      open(my $metricsFh, $file) or die "ERROR: could not read file $file: $!";
      my $header = <$metricsFh>;
      my $values = <$metricsFh>;
      chomp($header, $values);
      close $metricsFh;

      # Header and value match up
      my @header = split(/\t/, $header);
      my @value  = split(/\t/, $values);
      for(my $i=0;$i<@header;$i++){
        $metric{$header[$i]} = $value[$i];
      }
      # Remove the directory and suffix of the filename
      $metric{File} = basename($metric{File}, qw(.shovill.skesa.fasta .shovill.skesa.gbk .shovill.skesa.bam));
    }
    # Combine the metrics
    # Alphabetize the header but take care of Filename separately
    my $filename = $metric{File};
    my @header = sort keys(%metric);
    @header = grep{!/File/} @header;

    # Write combined metrics to file
    open(my $combinedFh, ">", $metricsOut) or die "ERROR: could not write to $metricsOut: $!";
    print $combinedFh join("\t", "File", @header)."\n";
    print $combinedFh $metric{File};
    for my $h(@header){
      print $combinedFh "\t" . $metric{$h};
    }
    print $combinedFh "\n";
    close $combinedFh;

    $numMetrics++;
  }
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
    command("shovill --outdir $outdir --R1 $R1 --R2 $R2 --ram 64 --assembler skesa --cpus $$settings{numcpus} --keepfiles");
  };
  if($@){
    logmsg "shovill failed!\n$@";
    return "";
  }

  return $outdir
}

sub skesaToGraph{
  my($sample,$sampleInfo,$contigs,$settings)=@_;

  my($R1,$R2) = @{ $$sampleInfo{fastq} };
  my $outdir="$$settings{tempdir}/$sample";
  mkdir($outdir);

  logmsg "Creating gfa file for $sample";
  command("gfa_connector --cores $$settings{numcpus} --reads $R1 $R2 --contigs $contigs --gfa $outdir/contigs.gfa --csv $outdir/contigs.csv --kmer 75 > $outdir/gfa_connecter.log 2>&1");
  command("echo === $outdir/gfa_connecter.log ===; echo ...; tail $outdir/gfa_connecter.log");

  return $outdir;

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

