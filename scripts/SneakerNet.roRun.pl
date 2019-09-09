#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

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
  GetOptions($settings,qw(help createsamplesheet tempdir=s numcpus=i outdir=s)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  mkdir($$settings{tempdir}) if(!-e $$settings{tempdir});
  $$settings{numcpus}||=1;
  $$settings{outdir}||="sneakernet.out";

  die usage() if($$settings{help} || !@ARGV);

  for my $dir(@ARGV){
    my $sneakernetDir = makeSneakernetDir($dir,$settings);
    saveSneakernetDir($sneakernetDir, $$settings{outdir});
  }
  rmdir($$settings{tempdir}) || logmsg "WARNING: might not have been able to remove tempdir $$settings{tempdir}";

  return 0;
}

sub makeSneakernetDir{
  my($dir,$settings)=@_;

  my $outdir="$$settings{tempdir}/runData";
  mkdir $outdir;

  my @fastq       = glob("$dir/Data/Intensities/BaseCalls/*.fastq.gz");
  my $snok        = "$dir/snok.txt";
  my $sampleSheet = "$dir/Data/Intensities/BaseCalls/SampleSheet.csv";
  my $sampleSheet2= "$dir/SampleSheet.csv"; # backup in case the first is not found
  my $sampleSheet3= (glob("$dir/*.csv"))[0]; # backup in case the first is not found
  my $config      =  "$dir/Data/Intensities/BaseCalls/config.xml";
  my @interop     =  glob("$dir/InterOp/*");
  my @xml         = ("$dir/CompletedJobInfo.xml",
                     "$dir/runParameters.xml",
                     "$dir/GenerateFASTQRunStatistics.xml",
                     "$dir/RunInfo.xml",
                    );
  
  if(! -e $sampleSheet){
    if(! -e $sampleSheet2){
      if($$settings{createsamplesheet}){
        # do nothing
      }
      # If no samplesheet2, is there _any_ csv file?
      elsif(-e $sampleSheet3){
        logmsg "Found $sampleSheet3 and so I will use it.";
        $sampleSheet2 = $sampleSheet3;
      } else {
        die "ERROR: could not find the samplesheet in either $sampleSheet or $sampleSheet2";
      }
    }
    $sampleSheet = $sampleSheet2;
  }
  if(!@interop){
    @interop  = glob("$dir/QC/InterOp/*");
    if(!@interop){
      logmsg "ERROR: no interop files were found in $dir";
    }
  }
  #for(@fastq, $config, @interop, @xml){
  #  if(!-e $_){
  #    die "ERROR: file does not exist: $_";
  #  }
  #}

  if($$settings{createsamplesheet}){
    $sampleSheet = createSampleSheet($dir, $outdir, $settings);
  }

  if(-e $snok){
    cp($snok,"$outdir/".basename($snok));
  } else {
    logmsg "snok.txt not found. I will not read from it.";
    # "touch" the snok file
    open(my $fh, ">>", "$outdir/".basename($snok)) or die "ERROR: could not touch $outdir/".basename($snok).": $!";
    close $fh;
  }
  if(!@fastq){
    # If there aren't any fastq files, try to see if they're in the main dir
    @fastq = glob("$dir/*.fastq.gz $dir/*/*.fastq.gz");
    if(!@fastq){
      logmsg "WARNING: no fastq files were found; attempting bcl2fastq";
      @fastq=bcl2fastq($dir,$settings);
      if(!@fastq){
        die "ERROR: could not find any fastq files in $dir";
      }
    }
  }

  for(@fastq, $sampleSheet, $config){
    my $to="$outdir/".basename($_);
    cp($_, $to);
  }
  mkdir "$outdir/QC";
  for(@xml){
    my $to="$outdir/QC/".basename($_);
    cp($_, $to);
  }
  mkdir "$outdir/QC/InterOp";
  for(@interop){
    my $to="$outdir/QC/InterOp/".basename($_);
    cp($_, $to);
  }

  return $outdir;
}

sub bcl2fastq{
  my($dir,$settings)=@_;
  my $fastqdir="$$settings{tempdir}/bcl2fastq";
  mkdir($fastqdir);

  #command("bcl2fastq --input-dir $dir/Data/Intensities/BaseCalls --runfolder-dir $dir --output-dir $fastqdir --processing-threads $$settings{numcpus} --demultiplexing-threads $$settings{numcpus} --barcode-mismatches 1 >&2");
  command("bcl2fastq --input-dir $dir/Data/Intensities/BaseCalls --runfolder-dir $dir --output-dir $fastqdir --processing-threads $$settings{numcpus} --barcode-mismatches 1 --ignore-missing-bcls >&2");

  my @fastq=glob("$$settings{tempdir}/bcl2fastq/*.fastq.gz");
  return @fastq;
}

sub saveSneakernetDir{
  my($tmpdir,$outdir,$settings)=@_;
  system("mv -v $tmpdir $outdir 1>&2");
  die if $?;
  #File::Copy::mv($tmpdir,$outdir) or die "ERROR: could not move $tmpdir to $outdir: $!";
  return 1;
}

sub cp{
  my($from,$to)=@_;
  if(-e $to && -s $to > 0){
    logmsg "Found $to. Not copying";
    return 1;
  }
  logmsg "cp $from to $to";
  my $return = link($from, $to) ||
    File::Copy::cp($from,$to) or warn "ERROR: could not copy $from to $to: $!";
  open(my $fh, ">>", $to) or die "ERROR: could not write to $to: $!";
  close $fh;
  return $return;
}

sub createSampleSheet{
  my($dir, $outdir, $settings) = @_;

  my $numSamples=0;
  my $samplesheet = "$outdir/SampleSheet.csv";
  if(-e $samplesheet){
    die "ERROR: was going to create a samplesheet but it already exists at $samplesheet";
  }
  open(my $fh, ">", $samplesheet) or die "ERROR: could not write to $samplesheet: $!";
  print $fh "[Data]\n";
  for my $demuxSamples(glob("$dir/*_*.csv")){ # should just be one file
    open(my $demuxFh, "<", $demuxSamples) or die "ERROR: could not read from $demuxSamples: $!";
    my $found_the_samples = 0;
    while(<$demuxFh>){
      s/^\s+|\s+$//g; # whitespace trim

      if(/,Sample,/){
        s/,Sample,/,SampleID,/;
      }
      if(/Flowcell/){
        next;
      }
      if(/^,+$/){
        next;
      }
      if(/^\s*$/){
        next;
      }
      print $fh $_."\n"; # add back in newline
      $numSamples++;
    }
    close $demuxFh;
  }
  close $fh;

  # If no samples were printed, then just try to make a
  # sample sheet using fastq file information.
  if($numSamples == 0){
    $samplesheet = "$outdir/samples.tsv";
    logmsg "demultiplex sample sheets were not found. Creating one directly from fastq filenames into $outdir/samples.tsv";
    createSampleSheetOutOfThinAir($samplesheet, $dir, $settings);
  }

  return $samplesheet;
}

# Using lots of help from Taylor Griswold for creating the sample sheet
# out of thin air.
sub createSampleSheetOutOfThinAir{
  my($outSamplesheet, $indir, $settings) =@_;

  # get a list of fastq files for the next couple of steps
  # Do not include the full path; just the basename.
  my @fastq = map{basename($_)} glob("$indir/*.fastq.gz");
  
  #awk -F'_' '{print $1}' <(ls *.fastq.gz) | sed 's/-*$//g' | uniq > wgs-ids_detailed.txt 
  my %sample;
  for my $fastq(@fastq){
    my($sample) = split(/_/, $fastq);
    $sample =~ s/\-*$//; # remove trailing dashes
    $sample{$sample}=1;
  }
  my @sample = sort {$a cmp $b} keys(%sample);  # sample list is unique and sorted

  # for i in $(cat wgs-ids_detailed.txt); do R1=`realpath $(ls *$i*_R1*.fastq.gz)`; R2=`realpath $(ls *$i*_R2*.fastq.gz)`; wgsid=`basename $R1 ".fastq.gz" | cut -d '_' -f 1`; echo -e $wgsid"\ttaxon=Salmonella;route=calcengine\t"$R1";"$R2; done > M1234-19-004.tsv 
  open(my $outFh, ">", $outSamplesheet) or die "ERROR: could not write to $outSamplesheet: $!";
  for my $name(@sample){
    my $R1 = (grep { $_=~/$name.*_R1/; } @fastq)[0];
    my $R2 = $R1;
    $R2 =~ s/_R1/_R2/;

    if(!-e "$indir/$R1" || !-e "$indir/$R2"){
      die "ERROR: could not find both $R1 and $R2 in $indir!";
    }

    print $outFh join("\t", $name, ".", "$R1;$R2")."\n";
  }
  close $outFh;

  return 1;
}


sub usage{
  "'SneakerNet-read-only': Parses an unaltered Illumina run and formats it
  into something usable for SneakerNet

  Usage: $0 illuminaDirectory [illuminaDirectory2...]
  
  --numcpus             1
  --outdir             ''
  --tempdir            ''
  --createsamplesheet     Also create a SampleSheet.csv or samples.tsv as fallback
  "
}
