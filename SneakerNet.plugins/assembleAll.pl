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

our $VERSION = "3.0";
our $CITATION= "Assembly plugin by Lee Katz. Uses SHOvill.";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help tempdir=s debug numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
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
      'gfa_connector' => 'gfa_connector --version 2>/dev/null | grep gfa_connector',
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
  assembleAll($dir,$settings);
  #logmsg "Metrics can be found in $metricsOut";

  #my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{version=>$VERSION,
    "Minimum contig length for assembly metrics" => "500bp",
  });

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/assemblyMetrics.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/assemblyMetrics_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Assembly metrics\"\n";
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

sub assembleAll{
  my($dir,$settings)=@_;
  
  # Find information about each genome
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    logmsg "ASSEMBLE SAMPLE $sample";

    my $outdir="$dir/SneakerNet/assemblies/$sample";
    my $outassembly="$outdir/$sample.shovill.skesa.fasta";
    my $outgfa     ="$outdir/$sample.shovill.skesa.gfa";
    my $outgbk="$outdir/$sample.shovill.skesa.gbk";
    my $outgff="$outdir/$sample.shovill.skesa.gff";

    # Run the assembly
    if(!-e $outassembly){
      my $assemblyDir = assembleSample($sample,$info,$settings);
      my $assembly    = "$assemblyDir/contigs.fa";
      my $gfa         = "$assemblyDir/contigs.gfa";
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
      cp($gfa,      $outgfa)      or die "ERROR copying $gfa => $outgfa\n  $!";

      #mkdir "$outdir/prodigal"; # just make this directory right away
    }
    
    if(!-e "$outdir/quast/report.html"){
      logmsg "Running quast on $outassembly > $outdir/quast";
      command("quast $outassembly --glimmer --output-dir $outdir/quast --threads $$settings{numcpus} --rna-finding");
    }
    else{
      logmsg "Found $outdir/quast/report.html. Not rerunning";
    }

  }
  
}



sub assembleSample{
  my($sample,$sampleInfo,$settings)=@_;

  my($R1,$R2) = @{ $$sampleInfo{fastq} };

  if(!$R1 || !-e $R1){
    logmsg "Could not find R1 for $sample ($R1). Skipping.";
    return "";
  }
  if(!$R2 || !-e $R2){
    logmsg "Could not find R2 for $sample ($R2). Skipping.";
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
    command("gfa_connector --cores $$settings{numcpus} --reads $R1 $R2 --use_paired_ends --contigs $outdir/contigs.fa > $outdir/contigs.gfa");
  };
  if($@){
    logmsg "shovill failed!\n$@";
    return "";
  }

  return $outdir
}

sub usage{
  print "Assemble all genomes
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
";
  exit(0);
}

