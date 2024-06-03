#!/usr/bin/env perl
# Hello World example

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Bio::SeqIO;
use Bio::FeatureIO::gff;

use FindBin qw/$RealBin/;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "2.0";
our $CITATION = "Genotyping SneakerNet plugin by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  my @exe = qw(kma);
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe       => \@exe,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];

  if(! -d "$dir/SneakerNet"){
    mkdir "$dir/SneakerNet";
  }
  if(! -d "$dir/SneakerNet/forEmail"){
    mkdir "$dir/SneakerNet/forEmail";
  }
  if(! -d "$dir/SneakerNet/genotype"){
    mkdir "$dir/SneakerNet/genotype";
  }

  my $table = genotypeSamples($dir, $settings);

  my $rawMultiQC = makeMultiQC($dir, $settings);

  my $kmaVersion = `kma -v`;
  chomp($kmaVersion);
  recordProperties($dir,{exe => \@exe, version=>$VERSION, table=>$table, kma_version=>$kmaVersion, mqc=>$rawMultiQC});

  logmsg "Table can be found in $table";

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/genotype.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/genotype_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Genotyping\"\n";
  print $outFh "#description: \"These are the raw results that feed into other interpretation plugins like genotyping for Escherichia.<br />$plugin v$VERSION $docLink $pluginLink\"\n";
  print $outFh "#anchor: '$anchor'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  my $header = <$fh>;
  print $outFh "hit\t".$header;

  my $resultCounter = 0;
  while(<$fh>){
    next if(/^#/);
    chomp;

    # need to make sample names unique
    my($sample, @F) = split /\t/;
    my $hit = $sample. '-' . ++$resultCounter;

    print $outFh join("\t", $hit, $sample, @F)."\n";
  }
  close $fh;

  return $outtable;
}

sub genotypeSamples{
  my($dir, $settings) = @_;

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");

    # See if there is a genotype database
    my $genotypeDatabases = [];
    if($$info{taxonRules}{genotype}){
      # If there is a database for this taxon, then
      # check to make sure it has been downloaded.
      $genotypeDatabases = checkGenotypeDatabases($$info{taxonRules}{genotype}, $settings);
    }
    else {
      logmsg "No databases listed for sample $sample under taxon $$info{taxon}";
      logmsg "Therefore, I will not genotype $sample.";
      next;
    }

    # Genotype with the current databases
    my $outdir = "$dir/SneakerNet/genotype/$sample";
    mkdir($outdir);
    # Genotype one database at a time
    my $numDatabases = @$genotypeDatabases;
    for(my $i=0; $i<$numDatabases; $i++){
      my $R1 = $$info{fastq}[0];
      my $R2 = $$info{fastq}[1];
      my $target = "$outdir/".basename($$genotypeDatabases[$i]);
      if(! -e "$target.res" || $$settings{force}){
        command("kma -ipe $R1 $R2 -o $target -t_db $$genotypeDatabases[$i] -ID 60.0");
      }
    }
  }

  # Make an output summary table
  my $table = "$dir/SneakerNet/forEmail/genotype.tsv";
  command("head -n 1 -q $dir/SneakerNet/genotype/*/*.res | sed 's/^#//' | head -n 1 | sed 's/^/sample\t/' > $table");
  for my $subdir(glob("$dir/SneakerNet/genotype/*")){
    next if(!-d $subdir);
    my $sample = basename($subdir);
    command("tail -n +2 -q $subdir/*.res | sort | uniq | sed 's/^/$sample\t/' >> $table");
  }
  return $table;
}
    
# Download any databases
sub checkGenotypeDatabases{
  my($urls, $settings) = @_;
  my @db;

  if(!-e "$RealBin/../db/genotype"){
    mkdir("$RealBin/../db/genotype");
  }

  $urls = (ref($urls) eq 'ARRAY')?$urls:[$urls];

  for my $u(@$urls){
    my $target = "$RealBin/../db/genotype/".basename($u);

    # Download if the fasta is not there yet
    if(!-e $target){
      command("wget '$u' -O $target");
    }

    # Index with kma
    if(!-e "$target.length.b"){
      command("kma index -i $target");
    }

    push(@db, $target);
  }

  return \@db;
}

sub usage{
  print "Run genotyping. Fasta databases can be listed in taxonProperties.conf.

  Usage: $0 MiSeq_run_dir
  --numcpus  1  Number of CPUs
  --force       Whether to overwrite results
  --debug       Print debugging information
  --version     Print the version and exit
  --help        This usage menu
  --tempdir  '' Specify the temporary directory
  ";
  exit 0;
}

