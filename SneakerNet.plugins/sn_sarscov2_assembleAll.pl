#!/usr/bin/env perl
# Assemble genomes for SARS-CoV-2 amplicons

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;
use List::Util qw/min max/;

use threads;
use Thread::Queue;

use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readTsv exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.1";
our $CITATION= "SARS-CoV-2 assembly plugin by Lee Katz.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      #'run_prediction_metrics.pl (CG-Pipeline)'     => "echo CG Pipeline version unknown",
      'run_assembly_metrics.pl (CG-Pipeline)'       => "echo CG Pipeline version unknown",
      cat                             => 'cat --version | head -n 1',
      wget                            => 'wget --version | head -n 1',
      samtools                        => 'samtools 2>&1 | grep Version:',
      bcftools                        => 'bcftools 2>&1 | grep Version:',
      cutadapt                        => 'cutadapt --version',
      seqtk                           => 'seqtk 2>&1 | grep -m 1 Version:',
      bgzip                           => 'bgzip --version | head -n1',
      tabix                           => 'tabix --version | head -n1',
      'v-annotate.pl (VADR)'          => 'v-annotate.pl -h | grep -m 1 [0-9]'
    }, $settings,
  );

  die usage() if($$settings{help} || !@ARGV);
  $$settings{VADRMODELDIR} ||= die "ERROR: VADRMODELDIR needs to be defined under config/settings.conf.
  Please see for more help: https://github.com/nawrockie/vadr/wiki/Coronavirus-annotation#how-to-annotate-sars-cov-2-sequences-with-vadr-1";
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

  recordProperties($dir,{version=>$VERSION,table=>$metricsOut});

  return 0;
}

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  my %sampleMetrics = ();
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    logmsg "ASSEMBLE SAMPLE $sample";

    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly="$outdir/$sample.bowtie2.bcftools.fasta";
    my $outgbk="$outdir/$sample.bowtie2.bcftools.gbk";

    # Run the assembly
    if(!-e $outassembly){
      my $assemblyDir = assembleSample($sample,$info,$settings);
      my $assembly    = "$assemblyDir/consensus.fasta";
      next if(! -d $assemblyDir);
      next if(! -e $assembly);
      next if(! -s $assembly);

      # Save the assembly
      mkdir $outdir;
      mkdir "$outdir/consensus";
      for my $srcFile(glob("$assemblyDir/*")){
        my $target = "$outdir/consensus/".basename($srcFile);
        cp($srcFile, $target) or die "ERROR copying $srcFile => $target\n  $!";
      }
      cp($assembly, $outassembly) or die "ERROR copying $assembly => $outassembly\n  $!";
    }

    # Assembly metrics
    my $asmMetricsFile = "$outdir/assemblyMetrics.tsv";
    if(! -e $asmMetricsFile){
      logmsg "Assembly metrics for $outassembly";
      command("run_assembly_metrics.pl --allMetrics --numcpus $$settings{numcpus} $outassembly >> $asmMetricsFile.tmp");
      mv("$asmMetricsFile.tmp", $asmMetricsFile) or die "ERROR: could not make $asmMetricsFile";
    }

    # Genome annotation
    #logmsg "PREDICT GENES FOR SAMPLE $sample";
    my $vadrPass = "$outdir/vadr/vadr.vadr.pass.tbl";
    my $vadrFail = "$outdir/vadr/vadr.vadr.fail.tbl";
    if(!-e $vadrPass || !-e $vadrFail){
      logmsg "ANNOTATE SAMPLE $sample";
      my $tmpVadrDir=annotateFasta($sample,$outassembly,$settings);
      mkdir("$outdir/vadr");
      for(glob("$tmpVadrDir/*")){
        next if(!-f $_);
        cp($_, "$outdir/vadr/".basename($_)) or die "ERROR copying $_ to $outdir/vadr/: $!";
      }
    }
    logmsg "Annotations can be found in $vadrPass and $vadrFail";

    my $cdsMetricsFile = "$outdir/predictionMetrics.tsv";
    if(!-e $cdsMetricsFile){
      logmsg "Prediction metrics for $vadrPass";
      my $passCds = readTbl($vadrPass, $settings);
      #print Dumper $passCds;
      logmsg "Prediction metrics for $vadrFail";
      my $failCds = readTbl($vadrFail, $settings);
      #print Dumper $failCds;

      my $numPassed = $$passCds{numFeatures}{mat_peptide} || 0;
      my $numFailed = $$failCds{numFeatures}{mat_peptide} || 0;

      open(my $fh, ">", $cdsMetricsFile.".tmp") or die "ERROR: could not write to $cdsMetricsFile.tmp: $!";
      print $fh join("\t", qw(File passCds failCds))."\n";
      print $fh join("\t", $outassembly, $numPassed, $numFailed)."\n";
      close $fh;
      mv("$cdsMetricsFile.tmp", $cdsMetricsFile) or die "ERROR: could not mv $cdsMetricsFile.tmp => $cdsMetricsFile: $!";
    }
    logmsg "Summary CDS info in $cdsMetricsFile";

    # Combine metrics
    my $annMetrics = readTsv($cdsMetricsFile, $settings);
    my $asmMetrics = readTsv($asmMetricsFile,$settings);
    my %asmAnnMetrics = (%$annMetrics, %$asmMetrics);
    while(my($filePath, $metricsHash) = each(%asmAnnMetrics)){
      while(my($metric, $value) = each(%$metricsHash)){
        $sampleMetrics{$sample}{$metric} = $value;
      }
      
      ## Add more metrics for this sample
      # percentage of Ns
      my %ntCounter;
      my $totalNt = 0;
      my $seqin = Bio::SeqIO->new(-file=>$outassembly);
      while(my $seq = $seqin->next_seq){
        my $sequence = $seq->seq;
        while($sequence =~ /(.)/g){
          $ntCounter{uc($1)}++;
        }
        $totalNt += length($sequence);
      }
      $seqin->close;

      $totalNt ||= ~0; # avoid divide by zero error by setting this number to something really high if zero.
      $sampleMetrics{$sample}{percentNs} = sprintf("%0.2f", $ntCounter{N} / $totalNt);
    }

  }

  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  if(!-e $metricsOut || $$settings{force}){
    logmsg "Combining metrics into $metricsOut";
    open(my $fh, ">", $metricsOut.".tmp") or die "ERROR: could not write to $metricsOut.tmp: $!";
    my @sample = sort {$a cmp $b} keys(%sampleMetrics);
    my @header = grep {!/File/} sort{$a cmp $b} keys($sampleMetrics{$sample[0]});
    print $fh join("\t", "File", @header) ."\n";
    for my $sample(@sample){
      print $fh basename($sampleMetrics{$sample}{File});
      for(@header){
        print $fh "\t".$sampleMetrics{$sample}{$_}
      }
      print $fh "\n";
    }
    close $fh;
    mv("$metricsOut.tmp", $metricsOut) or die "ERROR: could not move metrics file to $metricsOut: $!";
  }
  
  return $metricsOut;
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
  # TODO look at read metrics instead
  my $rawReadsFilesize = (stat($R1))[7] + (stat($R2))[7];
  if($rawReadsFilesize < 1000){
    if($$settings{force}){
      logmsg "Fastq files are tiny, but forcing because of --force.";
    } else {
      logmsg "Fastq files are too tiny. Skipping. Force with --force.";
      return "";
    }
  }

  my $outdir="$$settings{tempdir}/$sample";
  
  system("rm -rf '$outdir'"); # make sure any previous runs are gone
  mkdir $outdir;

  my $ref = $$sampleInfo{taxonRules}{reference_fasta};

  # Check for bowtie2 index files
  if(!-e "$ref.1.bt2"){
    my $log = "$ref.bowtie2-build.log";
    logmsg "Index not found. bowtie2-build log can be found at $log";
    command("bowtie2-build $ref $ref 2> $log");
  }

  ($R1, $R2) = adapterTrim($R1, $R2, $settings);

  logmsg "Mapping reads from $sample to $ref";
  my $samtoolsThreads = $$settings{numcpus}-1;
  my $sam       = "$outdir/unsorted.sam";
  my $sortedBam = "$outdir/sorted.bam";
  my $mpileup   = "$outdir/samtools.mpileup";
  command("bowtie2 --sensitive-local -p $$settings{numcpus} -x $ref -1 $R1 -2 $R2 -S $sam 2>&1");
  command("samtools view -b $sam | samtools sort --threads $samtoolsThreads - -o $sortedBam 2>&1");

  logmsg "Calling variants from mapped reads";
  my $vcf = "$outdir/out.vcf";
  system("samtools mpileup -aa -d 0 -uf $ref $sortedBam > $mpileup");
  die if $?;
  system("bcftools call -Mc < $mpileup > $vcf");
  die if $?;
  command("bgzip $vcf");
  command("tabix $vcf.gz");

  my $mindepth = $$sampleInfo{taxonRules}{coverage};
  my $consensusfasta = "$outdir/consensus.fasta";
  system("perl $RealBin/helper/vcf_mask_lowcoverage.pl --bam $sortedBam --vcf $vcf.gz --depth $mindepth --reference $ref --consout $consensusfasta");
  die if $?;

  # Some cleanup of needlessly large files that we just
  # want to make sure are gone
  unlink($sam);
  unlink("$outdir/out.vcf.gz");
  unlink("$outdir/out.vcf.gz.tbi");

  return $outdir
}

# annotate with vadr
sub annotateFasta{
  my($sample,$assembly,$settings)=@_;
  my $outdir="$$settings{tempdir}/$sample/vadr";
  system("rm -rf $outdir");

  mkdir "$$settings{tempdir}/$sample";

  command("v-annotate.pl -r -s --nomisc --lowsimterm 2 --mxsize 64000 --mdir $$settings{VADRMODELDIR} --mkey NC_045512 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf --lowsc 0.75 $assembly $outdir");

  return $outdir;
}

sub adapterTrim{
  my($R1in, $R2in, $settings) = @_;
  #my($read1, $read2, $settings) = @_;
  my $threads = $$settings{numcpus};

  #my $R1in  = "$$settings{tempdir}/in.".basename($read1);
  #my $R2in  = "$$settings{tempdir}/in.".basename($read2);
  my $R1out = "$$settings{tempdir}/".basename($R1in);
  my $R2out = "$$settings{tempdir}/".basename($R2in);

  #cp($read1, $R1in) or die $!;
  #cp($read2, $R2in) or die $!;

  my $log = "$$settings{tempdir}/cutadapt.log";
  open(my $logFh, ">", $log) or die "ERROR: could not write to $log: $!";
  close $logFh;

  logmsg "Cutadapt log file: $log";

  command("cutadapt -j $threads -g GTTTCCCAGTCACGATA -G GTTTCCCAGTCACGATA -a TATCGTGACTGGGAAAC -A TATCGTGACTGGGAAAC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -G ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -n 3 -m 75 -q 25 --interleaved $R1in $R2in 2>>$log | cutadapt -j $threads --interleaved -m 75 -u 30 -u -30 -U 30 -U -30 -o $R1out -p $R2out - >> $log 2>&1");

  return($R1out, $R2out);
}

# Read tbl format into a hash of features
# https://www.ncbi.nlm.nih.gov/projects/Sequin/table.html
sub readTbl{
  my($infile, $settings) = @_;

  # Variable for just feature information
  my %feature;
  # Variable for containing all features + meta information
  my %uber;
  # Contains counts of features by $featureKey
  my %numFeatures;

  my $currentSeqid = "UNKNOWN";
  my $currentTableName = "UNKNOWN";
  my ($featureStart, $featureStop, $featureKey);
  open(my $fh, "<", $infile) or die "ERROR: could not read $infile: $!";
  while(<$fh>){
    chomp;

    my @F = split(/\s+/, $_);
    next if(/^\s*$/);

    # get the contig information if it starts with >Feature
    if($F[0] =~ /^>Feature/){
      shift(@F); # Removes ">Feature" from array
      ($currentSeqid, $currentTableName) = @F;
      die "ERROR: it seems like there should be a sequence identifer but there wasn't one in $infile" if(!$currentSeqid);
      next;
    }
    
    # Get gene location information
    #   This is a three column format with start/stop/featureKey
    if(defined($F[0]) && $F[0] ne ""){
      
      # In VADR, there can be a line after all features saying
      #   Additional note(s) to submitter
      # This clues us into that all the rest of the text is not
      # part of the format.
      if($F[0] =~ /^\s*Additional note/i){
        last;
      }

      ($featureStart, $featureStop, $featureKey) = @F;
      if(!defined($featureKey)){
        $featureKey = "UNTYPED_FEATURE";
      }
      $numFeatures{$featureKey}++;
      next;
    }

    # Get qualifier key and value.
    # This line type is indented with threee empty columns.

    my($qualifierKey, $qualifierValue) = grep {/./} @F;

    # This is the structure of the features hash.
    # It's not necessary to explicitly define it like this but
    # it's also possibly confusing.
    $feature{$currentSeqid} //= {};
    $feature{$currentSeqid}{$featureStart}{$featureStop} //= {};
    $feature{$currentSeqid}{$featureStart}{$featureStop}{$featureKey} //= {};

    $feature{$currentSeqid}{$featureStart}{$featureStop}{$featureKey}{$qualifierKey} = $qualifierValue;

    #print Dumper [[$featureStart, $featureStop, $featureKey], [$qualifierKey, $qualifierValue]];
  }

  %uber = (
    features    => \%feature,
    numFeatures => \%numFeatures,
  );

  return \%uber;
}

sub usage{
  print "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
";
  exit(0);
}

