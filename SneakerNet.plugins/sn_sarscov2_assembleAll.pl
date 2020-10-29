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

use threads;
use Thread::Queue;

use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.0";
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
      #'v-annotate.pl (VADR)'          => 'v-annotate.pl -h | grep -m 1 [0-9]'
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

      mkdir "$outdir/vadr"; # just make this directory right away
    }

    # Genome annotation
    #logmsg "PREDICT GENES FOR SAMPLE $sample";
    #if(!-e $outgbk){
    #  my $gbk=annotateFasta($sample,$outassembly,$settings);
    #  cp($gbk,$outgbk) or die "ERROR: could not copy $gbk to $outgbk: $!";
    #}

  }
  
  # run assembly metrics with min contig size=0.5kb
  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  logmsg "Running metrics on the genbank files at $metricsOut";

  my @thr;
  my $Q=Thread::Queue->new(glob("$dir/SneakerNet/assemblies/*/*.bowtie2.bcftools.fasta"));
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&assemblyMetricsWorker,$Q,$settings);
    $Q->enqueue(undef);
  }
  for(@thr){
    $_->join;
  }
  
  command("cat $$settings{tempdir}/worker.*/metrics.tsv | head -n 1 > $metricsOut"); # header
  command("sort -k1,1 $$settings{tempdir}/worker.*/metrics.tsv | uniq -u >> $metricsOut"); # content
  
  return $metricsOut;
}

sub assemblyMetricsWorker{
  my($Q,$settings)=@_;
  my $tempdir=tempdir("worker.XXXXXX", DIR=>$$settings{tempdir}, CLEANUP=>1);
  my $assemblyOut   ="$tempdir/assemblyMetrics.tsv";
  my $metricsOut    ="$tempdir/metrics.tsv";  # Combined metrics
  my $coverageOut   ="$tempdir/effectiveCoverage.tsv";
  command("touch $assemblyOut");

  my $numMetrics=0;
  while(defined(my $fasta=$Q->dequeue)){
    # Bring the shovill bam file over
    my $bam = dirname($fasta)."/consensus/sorted.bam";
    
    logmsg "asm metrics for $fasta";
    command("run_assembly_metrics.pl --allMetrics --numcpus 1 $fasta >> $assemblyOut");

    logmsg "Effective coverage for $bam";
    command("samtools depth -aa $bam > $coverageOut");
    my $total=0;
    my $numSites=0;
    for my $line(`cat $coverageOut`){
      chomp($line);
      my(undef, undef, $depth) = split(/\t/, $line);
      $total+=$depth;
      $numSites++;
    }
    open(my $tmpFh, ">", $coverageOut) or die "ERROR: could not write to $coverageOut: $!";

    # Effective coverage: avoid divide by zero error.
    my $effectiveCoverage = "0.00";
    if($numSites > 0){
      $effectiveCoverage = sprintf("%0.2f", $total/$numSites);
    }

    print $tmpFh join("\t", qw(File effectiveCoverage))."\n";
    print $tmpFh join("\t", $fasta, $effectiveCoverage)."\n";
    close $tmpFh;
    $numMetrics++;
  }

  # Don't do any combining if no metrics were performed
  return if($numMetrics==0);

  # Combine the files
  open(my $assemblyFh, $assemblyOut) or die "ERROR: could not read $assemblyOut: $!";
  open(my $coverageFh, $coverageOut) or die "ERROR: could not read $coverageOut: $!";
  open(my $metricsFh, ">", $metricsOut) or die "ERROR: could not write to $metricsOut: $!";

  # Get and paste the header into the output metrics file
  my $assemblyHeader=<$assemblyFh>;
  my $coverageHeader=<$coverageFh>;
  chomp($assemblyHeader,$coverageHeader);
  my @assemblyHeader=split(/\t/,$assemblyHeader);
  my @coverageHeader=split(/\t/, $coverageHeader);
  print $metricsFh "File\tgenomeLength";
  for(@assemblyHeader, @coverageHeader){
    next if($_ =~ /File|genomeLength/);
    print $metricsFh "\t".$_;
  }
  # adding on percent Ns to the header
  print $metricsFh "\t"."percentNs";
  print $metricsFh "\n";

  # Get and paste the metrics from asm and pred into the output file
  while(my $assemblyLine=<$assemblyFh>){
    my $coverageLine=<$coverageFh>;
    chomp($assemblyLine,$coverageLine);

    my %F;
    @F{@assemblyHeader}=split(/\t/,$assemblyLine);
    @F{@coverageHeader}=split(/\t/,$coverageLine);
    # In case there was no assembly or predictions, add in a default value
    for(@assemblyHeader, @coverageHeader){
      $F{$_} //= 'NA';
    }
    # Also take care of genomeLength specifically
    $F{assembly} = $F{File}; # save the path to the asm
    $F{genomeLength} //= "NA";
    $F{File}=basename($F{assembly},qw(.gbk .fasta .bam));

    # Check percent Ns in the assembly
    my %ntCounter;
    my $totalNt = 0;
    my $seqin = Bio::SeqIO->new(-file=>$F{assembly});
    while(my $seq = $seqin->next_seq){
      my $sequence = $seq->seq;
      while($sequence =~ /(.)/g){
        $ntCounter{uc($1)}++;
      }
      $totalNt += length($sequence);
    }
    $seqin->close;
    $F{percentN} = sprintf("%0.2f", $ntCounter{N} / $totalNt);

    # Combine all the values into the metrics file
    print $metricsFh $F{File}."\t".$F{genomeLength};
    for(@assemblyHeader, @coverageHeader){
      next if($_ =~ /File|genomeLength/);
      print $metricsFh "\t".$F{$_}
    }
    # add percent Ns to the metrics
    print $metricsFh "\t".$F{percentN};

    print $metricsFh "\n";
  }
  close $_ for($assemblyFh,$coverageFh,$metricsFh);

  #system("ls $metricsOut; cat $metricsOut");
  #die "Problem with $metricsOut" if $?;
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

  # Some cleanup
  unlink($sam);
  unlink("$outdir/out.vcf.gz");
  unlink("$outdir/out.vcf.gz.tbi");

  return $outdir
}

# annotate with vadr
sub annotateFasta{
  my($sample,$assembly,$settings)=@_;

  # Ensure a clean slate
  my $outdir="$$settings{tempdir}/$sample/vadr";
  system("rm -rf $outdir");

  mkdir "$$settings{tempdir}/$sample";
  mkdir $outdir;

  my $outgff="$outdir/vadr.gff";
  my $outgbk="$outdir/vadr.gbk";

  logmsg "Annotating genes on $sample with VADR";
  eval{
    #command("prodigal -q -i $assembly -o $outgff -f gff -g 11 1>&2");
    # v-annotate.pl -s --nomisc --lowsimterm 2 --mxsize 64000 --mdir $VADRMODELDIR --mkey NC_045512 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf --lowsc 0.75 SRR12754173.bowtie2.bcftools.fasta vadr
    command("v-annotate.pl -s --nomisc --lowsimterm 2 --mxsize 64000 --mdir $$settings{VADRMODELDIR} --mkey NC_045512 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf --lowsc 0.75 $assembly $outdir/vadr");
  };
  # If there is an issue, push ahead with a zero byte file
  if($@){
    logmsg "There was an issue with predicting genes on sample $sample. This might be caused by a poor assembly.";
    open(my $fh, ">", $outgff) or die "ERROR: could not write to $outgff: $!";
    close $fh;
  }

  die;

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

sub usage{
  print "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
";
  exit(0);
}

