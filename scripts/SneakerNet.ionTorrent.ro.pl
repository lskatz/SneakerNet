#!/usr/bin/env perl
# Reads a read-only full version of an ion torrent run and
# turns it into a sneakernet folder.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/mv/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;
use JSON ();
use Encode qw/encode decode/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s numcpus=i outdir=s)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  $$settings{numcpus}||=1;
  $$settings{outdir}||="sneakernet.out";

  die usage() if($$settings{help} || !@ARGV);

  for my $dir(@ARGV){
    my $sneakernetDir = makeSneakernetDir($dir,$settings);
    saveSneakernetDir($sneakernetDir, $$settings{outdir});
  }
  rmdir($$settings{tempdir});

  return 0;
}

sub makeSneakernetDir{
  my($dir,$settings)=@_;

  my $outdir="$$settings{tempdir}/runData";
  mkdir $outdir;

  my $snok        = "$dir/snok.txt";
  my @remoteFastq = sort{(-s $a) <=> (-s $b) }
                      glob("$dir/plugin_out/FileExporter_out.*/*.fastq");

  my @localFastq;
  for my $old(@remoteFastq){
    my $new   = "$outdir/".basename($old);
    my $newGz = $new;
       $newGz =~ s/(\.gz)?$/.gz/;
    next if(! -f $old); # must be a file and not a symlink

    # Copy over the fastq file but make sure it is compressed
    if(! -e $newGz){
      if($old =~ /\.gz$/){
        cp($old, "$newGz.tmp");
      } else {
        logmsg "gzip $old => $newGz.tmp";
        system("gzip -c \Q$old\E > \Q$newGz.tmp\E");
        if($?){
          die "ERROR: could not gzip $old to $newGz.tmp";
        }
      }
      logmsg "  Removing .tmp from filename";
      mv("$newGz.tmp", $newGz)
        or die "ERROR renaming $newGz.tmp to $newGz - $!";
    }else{
      logmsg "SKIPPING: Found $newGz";
    }
    push(@localFastq, $newGz);
  }

  my $runInfo = iontorrentRunInfo($dir, $settings);
  my $barcodedSamples = $$runInfo{objects}[0]{barcodedSamples};
  my @sample = sort keys(%$barcodedSamples);
  my %barcode;
  for my $sample(@sample){
    while(my($barcode, $barcodeHash) = each(%{ $$barcodedSamples{$sample}{barcodeSampleInfo} })){
      $barcode{$barcode} = $sample;
    }
  }

  # Start off a sample sheet with columns: sampleName, [options], R1;R2
  # where R2 is not present if single end
  my $samplesheet = "$outdir/samples.tsv";
  open(my $fh, ">", $samplesheet) or die "ERROR writing to $samplesheet: $!";
  for my $fastq(@localFastq){
    my $longsample = basename($fastq, qw(.fastq.gz .fastq));
    my $sample = ".";
    my $barcode= ".";
    if($longsample =~ /([^\.]*).*(IonXpress_(\d+))/){
      $sample = $1;
      $barcode= $2;
    }
    # ok yes we have the sample name but in case it is listed
    # differently in the json file describing the run, use it
    # instead.
    if($barcode{$barcode}){
      $sample = $barcode{$barcode};
    }

    print $fh join("\t",
      $sample,
      "taxon=",
      $fastq,
    )."\n";
  }

  close $fh;

  if(-e $snok){
    cp($snok,"$outdir/".basename($snok));
  } else {
    logmsg "snok.txt not found. I will not read from it.";
    # "touch" the snok file
    open(my $fh, ">>", "$outdir/".basename($snok)) or die "ERROR: could not touch $outdir/".basename($snok).": $!";
    close $fh;
  }

  return $outdir;
}

sub iontorrentRunInfo{
  my($dir, $settings) = @_;

  my $jsonIn = "$dir/planned_run.json";
  if(! -e $jsonIn){
    die "ERROR: cannot find $jsonIn";
  }

  my $json=JSON->new;
  $json->utf8;           # If we only expect characters 0..255. Makes it fast.
  $json->allow_nonref;   # can convert a non-reference into its corresponding string
  $json->allow_blessed;  # encode method will not barf when it encounters a blessed reference
  $json->pretty;         # enables indent, space_before and space_after

  my $jsonStr = "";
  open(my $jsonFh, "<", $jsonIn) or die "ERROR reading $jsonIn: $!";
  {
    local $/ = undef;
    $jsonStr = <$jsonFh>;
  }
  close $jsonFh;

  # Need to check for valid utf8 or not
  eval{ my $strCopy = $jsonStr; decode('utf8', $strCopy, Encode::FB_CROAK) }
    or die "ERROR: $dir/planned_run.json yielded non-utf8 characters\nContents shown below:\n$jsonStr\n";

  my $runInfo = $json->decode($jsonStr);
  return $runInfo;
}

sub saveSneakernetDir{
  my($indir, $outdir)=@_;
  mkdir($outdir) if(! -e $outdir);
  for my $file(glob("$indir/*")){
    mv($file, "$outdir/")
      or logmsg "WARNING: could not move $file to $outdir/ - $!";
  }
}

sub cp{
  my($from,$to)=@_;
  if(-e $to && -s $to > 0){
    logmsg "Found $to. Not copying";
    return 1;
  }
  logmsg "cp $from to $to";
  my $return = link($from, $to) ||
    File::Copy::cp($from,$to) or die "ERROR: could not copy $from to $to: $!";
  open(my $fh, ">>", $to) or die "ERROR: could not write to $to: $!";
  close $fh;
  return $return;
}

sub usage{
  "Parses an unaltered Ion Torrent run and formats it
  into something usable for SneakerNet

  Usage: $0 iontorrentDirectory1 [iontorrentDirectory2...]
  
  --numcpus  1
  --outdir   ''
  "
}
