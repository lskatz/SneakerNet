#!/usr/bin/env perl
# Creates a pass/fail file for easy to read results

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp;
use FindBin;
use List::Util qw/sum/;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readTsv @rankOrder %rankOrder readKrakenDir exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "6.1";
our $CITATION="SneakerNet pass/fail by Lee Katz";

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version tempdir=s citation check-dependencies help debug force numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  if(!defined($$settings{kraken_contamination_min_rank})){
    logmsg "WARNING: kraken_contamination_min_rank was not set in settings.conf. Defaulting to G";
    logmsg "NOTE: Possible values for kraken_contamination_min_rank are single letters: @rankOrder";
    $$settings{kraken_contamination_min_rank} = 'G';
  }
  if(!defined($$settings{kraken_contamination_threshold})){
    logmsg "WARNING: kraken_contamination_threshold was not set in settings.conf. Defaulting to 25.0";
    $$settings{kraken_contamination_threshold} = 25.0;
  }

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";

  my $outfile=passfail($dir,$settings);
  logmsg "The pass/fail file is under $outfile";
  
  recordProperties($dir,{
      version=>$VERSION,table=>$outfile,
      kraken_contamination_min_rank  => $$settings{kraken_contamination_min_rank},
      kraken_contamination_threshold => $$settings{kraken_contamination_threshold},
  });

  return 0;
}

sub passfail{
  my($dir,$settings)=@_;

  my $failFile="$dir/SneakerNet/forEmail/passfail.tsv";
  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my $failHash=identifyBadRuns($dir,$sampleInfo,$settings);
  
  my @sample=keys(%$failHash);
  my @failHeader=sort keys(%{ $$failHash{$sample[0]} });

  open(my $failFh, ">", $failFile) or die "ERROR: could not write to $failFile: $!";
  print $failFh join("\t", "Sample", @failHeader)."\n";
  while(my($sample,$fail)=each(%$failHash)){
    print $failFh $sample;
    for(@failHeader){
      print $failFh "\t".$$fail{$_};
    }
    print $failFh "\n";
  }
  print $failFh "#1: fail\n#0: pass\n#-1: unknown\n";
  print $failFh "#Failure by Kraken is when the number of contaminant reads at rank $$settings{kraken_contamination_min_rank} or higher is >= $$settings{kraken_contamination_threshold}%\n";
  close $failFh;

  return $failFile;
}

sub identifyBadRuns{
  my($dir,$sampleInfo,$settings)=@_;

  my %whatFailed=();  # reasons why it failed

  # Read metrics
  my %readMetrics = ();
  # If the readMetrics.tsv file exists, read it.
  if(-e "$dir/readMetrics.tsv"){
    # Read the readMetrics file into %readMetrics
    open(READMETRICS,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
    my @header=split(/\t/,<READMETRICS>); chomp(@header);
    while(<READMETRICS>){
      chomp;
      my %F;
      @F{@header}=split(/\t/,$_);
      $F{File} = basename($F{File});
      $readMetrics{$F{File}} = \%F;
    }
    close READMETRICS;
  }
  # If the readmetrics file does not exist, then the values are blank
  else {
    logmsg "WARNING: readMetrics.tsv was not found. Unknown whether samples pass or fail by coverage and quality.";
  }

  # Assembly metrics
  my %assemblyMetrics = ();
  if(-e "$dir/SneakerNet/forEmail/assemblyMetrics.tsv"){
    my $assemblyMetricsTmp = readTsv("$dir/SneakerNet/forEmail/assemblyMetrics.tsv"); 

    # Remove extensions on the key to help with compatibility
    # down the line.
    # For example, either of these keys would be valid:
    #     SRR12530732.bowtie2.bcftools
    #     SRR12530732
    while(my($key, $metrics) = each(%$assemblyMetricsTmp)){
      my $newKey = $key;
      $newKey =~ s/\..*$//; # remove extension
      $assemblyMetrics{$key} = $metrics;
      $assemblyMetrics{$newKey} = $metrics;
    }
  }

  # Understand for each sample whether it passed or failed
  # on each category
  for my $samplename(keys(%$sampleInfo)){

    # Possible values for each:
    #  1: failed the category
    #  0: passed (did not fail)
    # -1: unknown
    my %fail=(
      coverage=>-1,
      quality => 0, # by default passes
      kraken  => -1,
      assembly=> -1,# collapse any issues with assembly into one category
    );

    # Skip anything that says undetermined.
    if($samplename=~/^Undetermined/i){
      next;
    }

    # Get the name of all files linked to this file through the sample.
    $$sampleInfo{$samplename}{fastq}//=[]; # Set {fastq} to an empty list if it does not exist
    my @file=@{$$sampleInfo{$samplename}{fastq}};

    # Get metrics from the fastq files
    my $totalCoverage = 0;
    my %is_passing_quality = (); # Whether a read passes quality
    for my $fastq(@file){
      my $fastqMetrics = $readMetrics{basename($fastq)};

      # Coverage
      $$fastqMetrics{coverage} //= '.';    # by default, dot.
      if($$fastqMetrics{coverage} eq '.'){ # dot means coverage is unknown
        $totalCoverage = -1; # -1 means 'unknown' coverage
      } else {
        $$fastqMetrics{coverage} ||= 0; # force it to be a number if it isn't already
        $totalCoverage += $$fastqMetrics{coverage};
        logmsg "Sample $samplename += $$fastqMetrics{coverage}x => ${totalCoverage}x" if($$settings{debug});
      }

      $$sampleInfo{$samplename}{taxonRules}{quality} ||= 0;
      # If avgQual is missing, then -1 for unknown pass status
      $$fastqMetrics{avgQuality} //= '.';
      if($$fastqMetrics{avgQuality} eq '.'){
        logmsg "$samplename/".basename($fastq)." qual is $$fastqMetrics{avgQuality}" if($$settings{debug});
        $is_passing_quality{$fastq} = -1;
      }
      # Set whether this fastq passes quality by the > comparison:
      # if yes, then bool=true, if less than, bool=false
      else {
        $is_passing_quality{$fastq} =  $$fastqMetrics{avgQuality} >= $$sampleInfo{$samplename}{taxonRules}{quality};
        $is_passing_quality{$fastq} += 0; # force to int/boolean
        logmsg "$samplename/".basename($fastq)." passes quality?  $is_passing_quality{$fastq} (boolean)" if($$settings{debug});
      }
    }

    # Set whether the sample fails coverage
    #logmsg "DEBUG"; $$sampleInfo{$samplename}{taxonRules}{coverage} = 5;
    if($totalCoverage < $$sampleInfo{$samplename}{taxonRules}{coverage}){
      logmsg "  I will fail this sample $samplename" if($$settings{debug});
      $fail{coverage} = 1;
    } else {
      logmsg "  I will not fail this sample $samplename" if($$settings{debug});
      $fail{coverage} = 0;
    }
    if($totalCoverage == -1){
      logmsg "  ==> -1" if($$settings{debug});
      $fail{coverage} = -1;
    }

    # if any filename fails quality, then overall quality fails
    while(my($fastq,$passed) = each(%is_passing_quality)){
      # If
      if($passed == 0){
        $fail{quality} = 1;
        logmsg "  I will fail sample $samplename on quality because ".basename($fastq) if($$settings{debug});
        last;
      }
      # Likewise if a file's quality is unknown, then it's all unknown
      if($passed == -1){
        $fail{quality} = -1;
        logmsg "  I will set quality failure as unknown for sample $samplename because ".basename($fastq) if($$settings{debug});
        last;
      }
    }

    # Did it pass by kraken / read classification
    # Start with a stance that we do not know (value: -1)
    $fail{kraken} = -1;
    my $krakenResults = readKrakenDir("$dir/SneakerNet/kraken/$samplename",{minpercent=>1});
    # If we have any results at all, then gain the stance
    # that it has not failed yet (value: 0)
    if(scalar(keys(%$krakenResults)) > 1){
      $fail{kraken} = 0;
    }
    for my $rank(reverse @rankOrder){
      # We can only call it conflicting if we are at the min rank or higher
      next if($rankOrder{$rank} < $rankOrder{$$settings{kraken_contamination_min_rank}});
      #next if(!$$krakenResults{$rank});

      # are there conflicting taxa at this rank?
      $$krakenResults{$rank} //= []; # Make sure this is defined before testing it in an array context
      my @taxon = @{$$krakenResults{$rank}};
      next if(scalar(@taxon) < 2);
      
      # Are the number of reads at least kraken_contamination_threshold?
      # The elements are sorted by percent and so the first element
      # is the dominant taxon.
      # Therefore we can judge contamination by the second element.
      # TODO instead, judge contamination by summing all percentages 
      #      except the zeroth.
      if($$krakenResults{$rank}[1]{percent} < $$settings{kraken_contamination_threshold}){
        next;
      }
      
      # If we made it through all filters, then we have reached
      # the point where we can call it contaminated.
      $fail{kraken} = 1;
      last; # no need to test any further ranks because it has failed
    }

    # Did it pass by kraken / read classification
    # Start with a stance that we do not know (value: -1)

    #die Dumper [\%assemblyMetrics, {%{$assemblyMetrics{$samplename}}}, $samplename];
    my $minNs = $$sampleInfo{$samplename}{taxonRules}{assembly_percentN};
    my $longestContigThreshold = $$sampleInfo{$samplename}{taxonRules}{longest_contig};
    if(defined($minNs) && defined($longestContigThreshold)){
      # Set up a subtest for the assembly.
      # If any of these assembly metrics fail, then
      # the assembly fails QC.
      # If all are unknown, then the assembly QC is unknown.
      # If all pass, then assembly QC passes.
      # Definitions:
      #   -1: unknown
      #    0: passes
      #    1: fails
      my %failAssembly = (
        minNs         => -1,
        longestContig => -1,
      );

      # QC for minimum number of Ns allowed
      my $percentNs = $assemblyMetrics{$samplename}{percentNs};
      if(!defined($percentNs)){
        $failAssembly{minNs} = -1;
      }
      elsif($percentNs > $minNs){
        $failAssembly{minNs} = 1;
      }
      elsif($percentNs <= $minNs){
        $failAssembly{minNs} = 0;
      }

      # QC for longest contig
      my $longestContig = $assemblyMetrics{$samplename}{longestContiguous};
      if(!defined($longestContig)){
        $failAssembly{longestContig} = -1;
      }
      elsif($longestContig >= $longestContigThreshold){
        $failAssembly{longestContig} = 0;
      }
      elsif($longestContig <  $longestContigThreshold){
        $failAssembly{longestContig} = 1;
      }

      # Compile assembly failure codes into a single code
      my $numFailed = 0;
      my $numTests  = 0;
      while(my($why, $failCode) = each(%failAssembly)){
        if($failCode==1){
          $numFailed++;
        }
        elsif($failCode==-1){
          next;
        }
        $numTests++;
      }
      if($numTests < 1){
        $fail{assembly} = -1;
      }
      else{
        if($numFailed > 0){
          $fail{assembly} = 1;
        } else {
          $fail{assembly} = 0;
        }
      }

    }
    # If assembly thresholds are not defined, then it is
    # unknown whether this sample passed.
    else {
      $fail{assembly} = -1;
    }

    # Save the culmination of pass/fail for this sample
    $whatFailed{$samplename} = \%fail;
  }


  return \%whatFailed;
}

sub usage{
  print "Passes or fails a run based on available information
  Usage: $0 MiSeq_run_dir
  --version
";
  exit(0);
}

