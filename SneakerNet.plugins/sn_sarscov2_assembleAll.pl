#!/usr/bin/env perl
# Assemble genomes for SARS-CoV-2 amplicons

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir tempfile/;
use File::Copy qw/mv cp/;
use File::Spec;
use Bio::SeqIO;
use Bio::FeatureIO::gff;
use List::Util qw/min max/;

use threads;
use Thread::Queue;

use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readTsv exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;
use GD; # for fonts
use GD::Graph::lines;

our $VERSION = "1.5";
our $CITATION= "SARS-CoV-2 assembly plugin by Lee Katz.";

# A message to show in the report if any
my $warningMsg="";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  my @exe = qw(run_assembly_metrics.pl cat wget samtools bcftools trimmomatic seqtk bgzip tabix bowtie2 bowtie2-build);
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

  my $imagePath = makeCoverageGraph($dir, $settings);

  recordProperties($dir,{
    version  => $VERSION,
    table    => $metricsOut,
    warnings => $warningMsg,
    image    => $imagePath,
    exe      => \@exe,
  });

  return 0;
}

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  my %sampleMetrics = ();
  my $cdsOverabundance = 0; # if any sample has > 100% CDS
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    logmsg "ASSEMBLE SAMPLE $sample";

    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly=File::Spec->rel2abs("$outdir/$sample.bowtie2.bcftools.fasta");
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

    my $cdsMetricsFile = "$outdir/predictionMetrics.tsv";
    if(!-e $cdsMetricsFile){
      my $tmppath = annotateFastaWithBioPerl($info, $outassembly, $settings);
      mv($tmppath, $cdsMetricsFile) or die $!;
    }
    logmsg "Summary CDS info in $cdsMetricsFile";

    # Combine metrics
    my $annMetrics = readTsv($cdsMetricsFile, $settings);
    my $asmMetrics = readTsv($asmMetricsFile,$settings);

    # Combine hashes into $$annMetrics
    while(my($filePath, $metricsHash) = each(%$asmMetrics)){
      while(my($metric, $value) = each(%$metricsHash)){
        $$annMetrics{$filePath}{$metric} = $value;
      }
    }

    while(my($filePath, $metricsHash) = each(%$annMetrics)){
      while(my($metric, $value) = each(%$metricsHash)){
        $sampleMetrics{$sample}{$metric} = $value;
      }
      
      ## Add more metrics for this sample
      my %ntCounter;
      my $totalNt = 0;
      my @contiguous;
      my $seqin = Bio::SeqIO->new(-file=>$outassembly);
      while(my $seq = $seqin->next_seq){
        my $sequence = $seq->seq;
        # percentage of Ns
        while($sequence =~ /(.)/g){
          $ntCounter{uc($1)}++;
        }
        
        # Longest contiguous section. GISAID(?) wants >10kb
        push(@contiguous, split(/[Nn]{1,}/, $sequence));
        $totalNt += length($sequence);
      }
      $seqin->close;
      my $longestContiguous = (sort{$b<=>$a} map{length($_)} @contiguous)[0];
      $sampleMetrics{$sample}{longestContiguous} = $longestContiguous;

      $totalNt ||= ~0; # avoid divide by zero error by setting this number to something really high if zero.
      $sampleMetrics{$sample}{percentNs} = sprintf("%0.2f", $ntCounter{N} / $totalNt);

      if($sampleMetrics{$sample}{expectedCdsPercentage} > 1){
        $cdsOverabundance = 1;
      }

      # remove any extraneous fields
      delete($sampleMetrics{$sample}{file});
    }

  }
  if($cdsOverabundance){
    $warningMsg .= "Some samples have reported >100% CDS. ";
  }

  my $metricsOut="$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  if(!-e $metricsOut || $$settings{force}){
    logmsg "Combining metrics into $metricsOut";
    open(my $fh, ">", $metricsOut.".tmp") or die "ERROR: could not write to $metricsOut.tmp: $!";
    my @sample = sort {$a cmp $b} keys(%sampleMetrics);
    my @header = grep {!/File/} sort{$a cmp $b} keys(%{ $sampleMetrics{$sample[0]} });
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

sub makeCoverageGraph{
  my($dir,$settings)=@_;

  my $png = "$dir/SneakerNet/forEmail/depth.png";
  my $maxPos = 0; # what is the max coordinate of any assembly
  my %depth; # depth{sample}[pos] = depthInt
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  my %sampleMetrics = ();
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");

    my $depthFile = "$dir/SneakerNet/assemblies/$sample/consensus/depth.tsv";

    # TODO account for genomes with multiple contigs
    # As of now, coordinates will restart with each new contig
    # which is not an issue if there is only one contig.
    my @depthOfContig = ();
    open(my $fh, $depthFile) or die "ERROR: could not read $depthFile: $!";
    while(<$fh>){
      chomp;
      my($contig, $pos, $depth) = split /\t/;
      $depthOfContig[$pos] = $depth;
    }
    close $fh;
    $maxPos = max($maxPos, scalar(@depthOfContig));
    
    $depth{$sample} = \@depthOfContig;
  }

  # Start a canvas for a line graph
  my $graphHeight = 300 + (10 * scalar(keys(%$sampleInfo)));
  logmsg "Graph height will be $graphHeight";
  my $graph = new GD::Graph::lines(1000, $graphHeight);
  $graph->set(
    title  => "Depth of reads per sample",
    x_label => "pos",
    y_label => "depth",
    correct_width => 0,
    x_tick_number => 'auto',
  );
  $graph->set_legend_font(GD::Font->Small);

  # get those samples into a set of arrays
  my @graphAmplitude;
  my @sampleName = sort keys %$sampleInfo;
  for(my $i=0;$i<@sampleName;$i++){
    $graphAmplitude[$i] = $depth{$sampleName[$i]};
  }

  $graph->set_legend(@sampleName);

  my $gd = $graph->plot([
    [1..$maxPos],
    @graphAmplitude,
  ]);
  if(!$gd){
    my $err = $graph->error;
    logmsg $err;
    $warningMsg .= "Could not graph depth into an image: $err ";
    return "";
  }
  
  open(my $pngFh, ">", $png) or die "ERROR: could not write to $png: $!";
  binmode($pngFh);
  print $pngFh $gd->png;
  close $pngFh;

  return $png;
}

sub assembleSample{
  my($sample,$sampleInfo,$settings)=@_;

  my($R1orig,$R2orig) = @{ $$sampleInfo{fastq} };

  if(!$R1orig || !-e $R1orig){
    logmsg "Could not find R1 for $sample. Skipping.";
    return "";
  }
  if(!$R2orig || !-e $R2orig){
    logmsg "Could not find R2 for $sample. Skipping.";
    return "";
  }

  # Get the reference ready
  my $ref = $$sampleInfo{taxonRules}{reference_fasta};
  # Check for bowtie2 index files
  if(!-e "$ref.1.bt2"){
    my $log = "$ref.bowtie2-build.log";
    logmsg "Index not found. bowtie2-build log can be found at $log";
    command("bowtie2-build $ref $ref 2> $log");
  }

  # Filter to reads that map to the reference
  my($R1,$R2) = readsThatMapTo($ref, $R1orig, $R2orig, $settings);

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

  logmsg "Trimming adapters";
  my $primersBed = $$sampleInfo{taxonRules}{primers_bed};
  if(!defined($primersBed)){
    die "ERROR: the primers_bed file is not defined. This is usually derived from the primers_bed_url key found in taxonProperties.conf.";
  }
  ($R1, $R2) = adapterTrim($R1, $R2, $primersBed, $ref, $settings);

  logmsg "Mapping reads from $sample to $ref";
  my $samtoolsThreads = $$settings{numcpus}-1;
  my $sam       = "$outdir/unsorted.sam";
  my $sortedBam = "$outdir/sorted.bam";
  my $mpileup   = "$outdir/samtools.mpileup";
  command("bowtie2 --sensitive-local -p $$settings{numcpus} -x $ref -1 $R1 -2 $R2 -S $sam 2>&1");
  command("samtools view --threads $samtoolsThreads -b $sam | samtools sort --threads $samtoolsThreads - -o $sortedBam 2>&1");
  command("samtools index $sortedBam");

  logmsg "Calling variants from mapped reads";
  my $vcf = "$outdir/out.vcf";
  system("samtools mpileup -aa -d 0 -uf $ref $sortedBam > $mpileup");
  die if $?;
  system("bcftools call -Mc < $mpileup > $vcf");
  die if $?;
  command("bgzip $vcf");
  command("tabix $vcf.gz");

  my $mindepth = $$sampleInfo{taxonRules}{per_base_coverage} || 1;
  my $consensusfasta = "$outdir/consensus.fasta";
  my $cmd = "perl $RealBin/helper/vcf_mask_lowcoverage.pl --bam $sortedBam --vcf $vcf.gz --depth $mindepth --reference $ref --consout $consensusfasta";
  logmsg "Using helper script\n  $cmd";
  system($cmd); die if $?;

  # Some cleanup of needlessly large files that we just
  # want to make sure are gone
  unlink($sam);
  unlink("$outdir/out.vcf.gz");
  unlink("$outdir/out.vcf.gz.tbi");

  command("samtools depth -aa $sortedBam > $outdir/depth.tsv");
  return $outdir;

}

# Filter to just reads that map to a reference
# The magic behind this is mapping and then samtools -F 4 -F 8
#   => read not unmapped and mate not unmapped
sub readsThatMapTo{
  my($ref, $R1orig, $R2orig, $settings) = @_;

  my $workingDir = File::Temp::tempdir("readsThatMapTo.XXXXXX", DIR=>$$settings{tempdir});

  my $samtoolsThreads = $$settings{numcpus}-1;
  my $sam = "$workingDir/bowtie2.sam";
  my ($fhR1, $R1new) = tempfile("R1.XXXXXX", SUFFIX=>".fastq", DIR=>$workingDir);
  my ($fhR2, $R2new) = tempfile("R2.XXXXXX", SUFFIX=>".fastq", DIR=>$workingDir);

  command("bowtie2 --very-fast-local -p $$settings{numcpus} -x $ref -1 $R1orig -2 $R2orig -S $sam 2>&1");
  command("samtools fastq -F 12 $sam -1 $R1new -2 $R2new --threads $samtoolsThreads");
  command("gzip $R1new");
  command("gzip $R2new");

  # Moderate cleanup before tempdir is destroyed anyway
  unlink($sam);

  return("$R1new.gz", "$R2new.gz");
}

sub annotateFastaWithBioPerl{
  my($info, $outassembly, $settings) = @_;

  my $tempdir = File::Temp::tempdir("annotateWithBioperl.XXXXXX",DIR=>$$settings{tempdir},CLEANUP=>1);

  # annotation strategy: load up the reference genome
  # with all its features. Replace the reference
  # genome sequence with the new sequence. See how the
  # CDS translations come out.
  # Early stop genes are pseudogenes.
  my %seq;
  my $refSeqCdsCounter = 0;
  my $refSeq = Bio::SeqIO->new(-file=>$$info{taxonRules}{reference_gbk});
  while(my $seq = $refSeq->next_seq){
    $seq{$seq->id} = $seq;

    for my $feat($seq->get_SeqFeatures){
      # Filter for only CDS
      next if($feat->primary_tag ne 'CDS');

      my $refProtObj = $feat->seq->translate;
      my $refAA = $refProtObj->seq;
      $refAA =~ s/\**$//; # remove stop codons at the end
      # If we don't see any stop codons, then let's
      # call it a CDS. A really naive method but it
      # should do the job for quick Q/C.
      # For example, the SARS-CoV-2 ref1ab gene will
      # probably not get counted here since it is a
      # weird virus gene that gets split after the
      # mRNA step.
      if($refAA !~ /\*/){
        $refSeqCdsCounter++;
      }
    }
  }
  $refSeq->close;

  my $altSeqIn = Bio::SeqIO->new(-file=>$outassembly);
  my $altSeqCdsCounter = 0;
  while(my $altSeq = $altSeqIn->next_seq){
    # replace the reference sequence with the alt
    # sequence but keep the features.
    my $seqid = $altSeq->id;
    # strip the NCBI version number
    $seqid =~ s/\.\d+$//;
    my $altSequence = $altSeq->seq;
    # add back in the reference features
    for my $feat($seq{$seqid}->get_SeqFeatures){
      $altSeq->add_SeqFeature($feat);
    }

    for my $feat($altSeq->get_SeqFeatures){
      # Filter for only CDS
      next if($feat->primary_tag ne 'CDS');

      my $altProtObj = $feat->seq->translate;
      my $altAA = $altProtObj->seq;
      $altAA =~ s/\**$//; # remove stop codons at the end
      if($altAA !~ /\*/){
        $altSeqCdsCounter++;
      }
    }

  }
  $altSeqIn->close;
  
  # Calculate percentage of CDSs intact.
  # Avoid divide by zero error by setting refCdsCount
  # to a really high number if not set.
  $refSeqCdsCounter ||= ~0;
  my $percentCds = sprintf("%0.2f", $altSeqCdsCounter/$refSeqCdsCounter);

  my $cdsMetricsFile = "$tempdir/prediction.tsv";
  open(my $fh, ">", $cdsMetricsFile) or die "ERROR writing to $cdsMetricsFile: $!";
  print $fh join("\t", qw(file refCdsCount altCdsCount expectedCdsPercentage))."\n";
  print $fh join("\t", $outassembly, $refSeqCdsCounter, $altSeqCdsCounter, $percentCds)."\n";
  close $fh;

  return $cdsMetricsFile;
}

# Pure perl adapter trimming with a bed file so that
# I don't have to use trimmomatic
sub adapterTrim{
  my($R1in, $R2in, $primersBed, $ref, $settings) = @_;
  my $threads = $$settings{numcpus};

  ## Need to get sequences from reference with bed coordinates
  # get reference contigs
  my $seqin = Bio::SeqIO->new(-file=>$ref);
  my %seq;
  while(my $seq=$seqin->next_seq){
    $seq{$seq->id} = $seq;
  }
  $seqin->close;

  # get correct sequences of primers
  my $PRIMERS = "$$settings{tempdir}/primers.fasta";
  open(my $seqout, '>', $PRIMERS) or die "ERROR: could not write to $PRIMERS: $!";
  open(my $bedFh, $primersBed) or die "ERROR: could not read from $primersBed: $!";
  #my @primerSeq;
  while(<$bedFh>){
    chomp;
    my($chrom, $start, $stop, $name, $score, $strand) = split(/\t/);
    my $seq = $seq{$chrom}->subseq($start+1,$stop);
    
    if($strand eq '-'){
      $seq = reverse($seq);
      $seq =~ tr/ACGTacgt/TGCAtgca/;
    }

    print $seqout ">$name\n$seq\n";
    #push(@primerSeq, $seq);
  }
  close($seqout);
  close($bedFh);

  #my $regexStr = join("|",@primerSeq);
  #my $regex = qr/^($regexStr)/i;

  my $R1out = "$$settings{tempdir}/".basename($R1in);
  my $R2out = "$$settings{tempdir}/".basename($R2in);

  my $MIN_BQ=3;
  my $TRIMOPT = "ILLUMINACLIP:$PRIMERS:1:30:11 LEADING:$MIN_BQ TRAILING:$MIN_BQ MINLEN:30 TOPHRED33";
  command("trimmomatic PE -threads $threads -phred33 \Q$R1in\E \Q$R2in\E \Q$R1out\E /dev/null \Q$R2out\E /dev/null $TRIMOPT 2>&1");

  return($R1out, $R2out);
}

## Trim a fastq file with just a compiled regex like
# $regex = qr/^(AAA|AAT|ATG)/
#sub trimOneFastq{
#  my($in, $out, $regex, $settings) = @_;
#
#  my $lineCounter = 0;
#  open(my $fh, "zcat $in |") or die "ERROR: could not open $in: $!";
#  open(my $fhOut, " | gzip -c > $out") or die "ERROR: could not pipe to $out: $!";
#  while(<$fh>){
#    $lineCounter++;
#    if($lineCounter % 4 == 2){
#      s/$regex//;
#    }
#    print $fhOut $_;
#  }
#  close $fhOut;
#  close $fh;
#  return 1;
#}
    

sub usage{
  print "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
";
  exit(0);
}

