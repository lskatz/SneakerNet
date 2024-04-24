#!/usr/bin/env perl
# Runs the NCBI "base-balance" metric which checks for whether
# there is approximately a 1 to 1 balance of As to Ts in the
# read set and a 1 to 1 balance between C and G.
# An offset of 1.3 is usually an indicator of problems.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/mv cp/;
use FindBin;

use threads;
use Thread::Queue;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig logmsg samplesheetInfo_tsv command/;
use List::MoreUtils qw/uniq/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

our $VERSION = "3.0";
our $CITATION = "Base Balance by Lee Katz";

# Global
my @nt=qw(A T C G N);
my @sampleHeader=("File",@nt,"A/T","C/G");

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies tempdir=s help inbox=s debug test force numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe => [ ],
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  
  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";

  my $out=baseBalanceAll($dir,$settings);
  
  my $target = "$dir/SneakerNet/forEmail/".basename($out);
  cp($out, $target) or die "ERROR: could not copy $out => $target";

  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{version=>$VERSION, table=>"$dir/SneakerNet/forEmail/".basename($out)});
  
  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/basebalance.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/bb_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Base balances\"\n";
  print $outFh "#description: Each A/T and C/G should equal about 1 and almost always between 0.7 and 1.3.<br />$plugin v$VERSION $docLink $pluginLink\n";
  print $outFh "#anchor: '$anchor'\n";
  print $outFh "#plot_type: 'bargraph'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  my $header = <$fh>;
  chomp($header);
  my @header = split(/\t/, $header);
  my %sample;
  while(<$fh>){
    next if(/^#/);
    my @F = split /\t/;
    my %F;
    @F{@header} = @F;
    $sample{$F{File}} = \%F;
  }
  close $fh;

  print $outFh join("\t", qw(sample ratio))."\n";
  while(my($sample, $values) = each(%sample)){
    print $outFh join("\t", $sample."-A/T", $$values{"A/T"})."\n";
    print $outFh join("\t", $sample."-C/G", $$values{"C/G"})."\n";
  }
  close $outFh;

  logmsg "Wrote $outtable";

  return $outtable;
}

sub baseBalanceAll{
  my($dir,$settings)=@_;

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my $outdir="$dir/SneakerNet/basebalance";

  mkdir $outdir; die if $?;
  my @outfile;
  my $outfile="$outdir/basebalance.tsv";

  my @thr;
  my $Q=Thread::Queue->new();
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&baseBalanceWorker,$outdir,$Q,$settings);
  }
  while(my($sampleName,$s)=each(%$sampleInfo)){
    next if(ref($s) ne 'HASH'); # avoid file=>name aliases

    $Q->enqueue($s);
  }
  for(@thr){
    $Q->enqueue(undef);
  }
  for(@thr){
    my $balanceFile=$_->join();
    push(@outfile,@$balanceFile);
  }
  
  open(OUTFILE,">",$outfile) or die "ERROR: could not open $outfile for writing: $!";
  print OUTFILE join("\t",@sampleHeader)."\n";
  close OUTFILE;
  command("grep -hv '^#' @outfile >> $outfile");

  return $outfile;
}

sub baseBalanceWorker{
  my($outdir,$Q,$settings)=@_;
  
  my @bbFile=();
  while(defined(my $sHash=$Q->dequeue)){
    my $bbFile=baseBalance($sHash,$outdir,$settings);
    push(@bbFile,$bbFile);
  }
  return \@bbFile;
}

sub baseBalance{
  my($sHash,$outdir,$settings)=@_;

  my $outfile="$outdir/$$sHash{sample_id}.tsv";;
  logmsg "Running base balance for $$sHash{sample_id} => $outfile";
  if(!$$settings{force}){
    return $outfile if(-e $outfile && (-s $outfile > 0));
  }

  open(my $outFh,">","$outfile.tmp") or die "ERROR: could not open $outfile.tmp for writing: $!";

  for my $fastq(@{ $$sHash{fastq} }){
    my $i=0;
    my %nt;
    $nt{$_}=1 for(@nt); # give a pseudocount just in case anything ends up being zero.
    open(my $fastqFh,"<",$fastq) or die "ERROR: could not open $fastq for reading: $!";
    while(<$fastqFh>){
      my $mod = ++$i % 4;
      if($mod==2){
        $_=uc($_);
        for my $nt(split(//)){
          $nt{$nt}++;
        }
      }
    }
    close $fastqFh;

    # header
    my $out='#'.join("\t",@sampleHeader)."\n";
    # values
    $out.=basename($fastq)."\t";
    $out.=$nt{$_}."\t" for(@nt);
    $out.=sprintf("%0.2f",($nt{A}/$nt{T}))."\t".sprintf("%0.2f",($nt{C}/$nt{G}))."\n";
    print $outFh $out;
  }
  print $outFh "# A/T is the number of As divided by Ts\n";
  print $outFh "# C/G is the number of Cs divided by Gs\n";
  print $outFh "# An expected ratio of randomly called nucleotides would be close to 1\n";
  close $outFh;

  mv("$outfile.tmp",$outfile) or die "ERROR: could not move $outfile.tmp to $outfile";
   
  return $outfile;
}


sub usage{
  print "Runs the base-balance metric
  Usage: $0 runDir
  --numcpus 1
  --version
  ";
  exit(0);
}

