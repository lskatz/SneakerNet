#!/usr/bin/env perl

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

our $VERSION = "1.0";
our $CITATION = "Genotyping SneakerNet plugin by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
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

  recordProperties($dir,{version=>$VERSION, table=>$table});

  logmsg "Table can be found in $table";

  return 0;
}

sub genotypeSamples{
  my($dir, $settings) = @_;

  my @subTable;
  
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);
  while(my($sample,$info)=each(%$sampleInfo)){
    next if(ref($info) ne "HASH");
    if($$info{taxon} ne 'Escherichia'){
      logmsg "SKIP: $sample taxon is $$info{taxon} and not Escherichia";
      next;
    }

    # To specifically genotype Escherichia, we need an O and H result
    # from sn_genotype.pl.
    my $oResFile = "$dir/SneakerNet/genotype/$sample/O_type.fsa.res";
    my $hResFile = "$dir/SneakerNet/genotype/$sample/H_type.fsa.res";
    if(! -e $oResFile){
      logmsg "SKIP: Escherichia sample $sample is missing a needed file $oResFile";
      next;
    }
    if(! -e $hResFile){
      logmsg "SKIP: Escherichia sample $sample is missing a needed file $hResFile";
      next;
    }
    my $oInfo = parseRes($oResFile, $settings);
    my $hInfo = parseRes($hResFile, $settings);

    my $O = (split(/_/, $$oInfo{template}))[-1];
    my $H = (split(/_/, $$hInfo{template}))[-1];

    my $table = "$dir/SneakerNet/genotype/$sample/genotype.Escherichia.tsv";
    logmsg "Writing to $table";
    open(my $fh, ">", $table) or die "ERROR: could not write to $table: $!";
    print $fh join("\t", qw(sample O H))."\n";
    print $fh join("\t", $sample, $O, $H)."\n";
    close $fh;

    push(@subTable, $table);
  }

  my $table = "$dir/SneakerNet/forEmail/genotype.Escherichia.tsv";
  open(my $fh, ">", $table) or die "ERROR: could not write to $table: $!";
  print $fh join("\t", qw(sample O H))."\n";
  for my $t(@subTable){
    open(my $subFh, $t) or die "ERROR: could not read sample table $t: $!";
    <$subFh>; # Burn the header
    my $line = <$subFh>;
    close $subFh;

    print $fh $line;
  }
  close $fh;
  
  return $table;
}

sub parseRes{
  my($resFile, $setting) = @_;
  
  my %bestResult;
  open(my $fh, $resFile) or die "ERROR: could not read $resFile: $!";
  my $header = <$fh>;
  $header =~ s/^#\s*//;
  chomp($header);
  $header = lc($header);
  my @header = split(/\t/, $header);
  while(<$fh>){
    chomp;
    my %F;
    my @F = split(/\t/, $_);
    for(my $i=0;$i<@header;$i++){
      $F{$header[$i]} = $F[$i];
    }

    if(!defined($bestResult{template})){
      %bestResult = %F;
    }

    # Determine if we need to update the best result
    if($F{score} > $bestResult{score}){
      %bestResult = %F;
    }
  }

  return \%bestResult;
}
    
sub usage{
  print "Parse genotyping specifically for Escherichia

  Usage: $0 run_dir
  --numcpus  1  Number of CPUs
  --force       Whether to overwrite results
  --debug       Print debugging information
  --version     Print the version and exit
  --help        This usage menu
  --tempdir  '' Specify the temporary directory
  ";
  exit 0;
}

