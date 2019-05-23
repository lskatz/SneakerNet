#!/usr/bin/env perl
# Reads a read-only full version of an ion torrent run and
# turns it into a sneakernet folder.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use FindBin;

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

  return 0;
}

sub makeSneakernetDir{
  my($dir,$settings)=@_;

  my $outdir="$$settings{tempdir}/runData";
  mkdir $outdir;

  my $snok        = "$dir/snok.txt";
  my @bam         = glob("$dir/exportedReports/*/basecaller_results/rawtf.basecaller.bam");
  #                                            ^-< this is the sample name

  # Start off a sample sheet with columns: sampleName, R1, R2
  # where R2 is not present if single end
  my $samplesheet = "$outdir/SampleSheet.tsv";
  open(my $fh, ">", $samplesheet) or die "ERROR writing to $samplesheet: $!";
  # TODO accept the three-column format in SneakerNet

  # Generate fastq files
  my $samtoolThreads = $$settings{numcpus} - 1;
  for my $bam(@bam){
    # Get some filenames
    my $sampleName;
    if($bam =~ m|^.+/exportedReports/(.+)/basecaller_results/|){
      $sampleName = $1;
    } else {
      die "ERROR: could not parse the sample name for $bam";
    }
    my $basename   = basename($bam, qw(.bam));
    my $fastq      = "$outdir/$sampleName.fastq";
    my $R2         = "."; # no idea how to define R2 in an ion torrent run

    # Fastq generation
    logmsg "Creating fastq from $bam";
    system("samtools fastq -@ $samtoolThreads -0 $fastq -N $bam");
    die "ERROR with $bam => $fastq" if $?;

    # Compress
    system("gzip -v $fastq"); die if $?;
    $fastq="$fastq.gz";

    # Samplesheet
    print $fh join("\t", $sampleName, basename($fastq),$R3)."\n";
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

sub saveSneakernetDir{
  my($indir, $outdir)=@_;
  system("mv -v '$indir' '$outdir'");
  die if $?;
}

sub cp{
  my($from,$to)=@_;
  if(-e $to && -s $to > 0){
    logmsg "Found $to. Not copying";
    return 1;
  }
  logmsg "cp $from to $to";
  my $return=File::Copy::cp($from,$to) or die "ERROR: could not copy $from to $to: $!";
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
