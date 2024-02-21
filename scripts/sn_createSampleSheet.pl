#!/usr/bin/env perl
# Creates a sample sheet for a folder

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Copy qw/mv/;
use File::Path qw/remove_tree/;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Find qw/find/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help fasta tempdir=s numcpus=i )) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  $$settings{numcpus}||=1;

  die usage() if($$settings{help} || !@ARGV);

  my $samplesheet = "$ARGV[0]/samples.tsv";
  createSampleSheetOutOfThinAir($samplesheet, $ARGV[0], $settings);

  return 0;
}

# Using lots of help from Taylor Griswold for creating the sample sheet
# out of thin air.
sub createSampleSheetOutOfThinAir{
  my($outSamplesheet, $indir, $settings) =@_;

  # Set the extension to look for as fastq.gz,
  # unless --fasta is set
  my $ext = "fastq.gz";
  if($$settings{fasta}){
    $ext="fasta";
  }

  # get a list of fastq files for the next couple of steps
  # Do not include the full path; just the basename.
  my %fastq; # basename => full path
  find({wanted=>sub{
    my $path = $File::Find::name;
    # Only look for files and not directories or symlinks
    return if(!-f $path);
    # Keep only fasta or fastq.gz files
    return if($path !~ /\.$ext$/);

    my $basename = basename($path);

    $fastq{$basename} = $path;
  }}, $indir);
  my @fastq = sort {$a cmp $b} keys(%fastq);
  
  #awk -F'_' '{print $1}' <(ls *.fastq.gz) | sed 's/-*$//g' | uniq > wgs-ids_detailed.txt 
  my %sample;
  if($$settings{fasta}){
    for my $fastq(@fastq){
      $sample{basename($fastq, ".fasta")}=1;
    }
  }
  else {
    for my $fastq(@fastq){
      my($sample) = split(/_/, $fastq);
      $sample =~ s/\-*$//; # remove trailing dashes
      $sample{$sample}=1;
    }
  }
  my @sample = sort {$a cmp $b} keys(%sample);  # sample list is unique and sorted

  # for i in $(cat wgs-ids_detailed.txt); do R1=`realpath $(ls *$i*_R1*.fastq.gz)`; R2=`realpath $(ls *$i*_R2*.fastq.gz)`; wgsid=`basename $R1 ".fastq.gz" | cut -d '_' -f 1`; echo -e $wgsid"\ttaxon=Salmonella;route=calcengine\t"$R1";"$R2; done > M1234-19-004.tsv 
  open(my $outFh, ">", $outSamplesheet) or die "ERROR: could not write to $outSamplesheet: $!";
  for my $name(@sample){
    if($$settings{fasta}){
      #my $taxon = estimateSpecies($fastq[0], $settings);
      print $outFh join("\t", $name, ".", $fastq[0])."\n";
    }
    else{
      my $R1 = (grep { $_=~/$name.*_R1|$name.*_1\.f.*q/; } @fastq)[0];
      my $R2 = $R1;
      $R2 =~ s/_R1/_R2/;

      # If the substitution didn't work, then it might be a _1.fastq.gz format
      if($R1 eq $R2){
        $R2 =~ s/_1\./_2\./;
      }

      # If R1 and R2 are still the same then we have problems
      if($R1 eq $R2){
        die "ERROR: could not find the pair for R1 $R1";
      }

      # Check that the keys exist
      if(!$fastq{$R1} && !$fastq{$R2}){
        die "ERROR: R1 and R2 are not both represented as $R1 and $R2";
      }
      # Check that the fastq files exist
      if(!-e $fastq{$R1}){
        die "ERROR: could not find R1 for $name at $fastq{$R1}";
      }
      if(!-e $fastq{$R2}){
        die "ERROR: could not find R2 for $name at $fastq{$R2}";
      }

      #if(!-e "$indir/$R1" || !-e "$indir/$R2"){
      #  die "ERROR: could not find both $R1 and $R2 in $indir!";
      #}
      print $outFh join("\t", $name, ".", "$R1;$R2")."\n";
    }

  }
  close $outFh;

  return 1;
}

#sub estimateSpecies{
#  my($seq, $settings) = @_;
#
#  system("mash dist $seq ~/bin/SneakerNet/db/fasta/*.fasta 2>/dev/null | sort -k3,3n");
#  die;
#}

sub usage{
  "Creates a samples.tsv out of thin air if there are fastq files

  Usage: $0 illuminaDirectory
  --fasta  Instead of fastq files, look for fasta assembly files
  "
}

