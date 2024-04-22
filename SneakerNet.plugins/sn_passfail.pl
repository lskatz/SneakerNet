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

our $VERSION = "8.0";
our $CITATION="SneakerNet pass/fail by Lee Katz";

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version verbose tempdir=s citation check-dependencies help debug force numcpus=i)) or die $!;
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
  
  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{
      version=>$VERSION,table=>$outfile,
      kraken_contamination_min_rank  => $$settings{kraken_contamination_min_rank},
      kraken_contamination_threshold => $$settings{kraken_contamination_threshold},
      mqc => $rawMultiQC,
  });

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/passfail.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/passfail_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"Pass/fail information\"\n";
  print $outFh "#description: \"$plugin v$VERSION $docLink $pluginLink\"\n";
  print $outFh "#anchor: '$anchor'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  while(<$fh>){
    next if(/^#/);
    chomp;

    my($sample, @F) = split /\t/;

    # Change integers to pass/fail strings
    for(@F){
      s/-1/unknown/;
      s/1/fail/;
      s/0/pass/;
    }
    print $outFh join("\t", $sample, @F)."\n";
  }
  close $fh;

  return $outtable;
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
      #$F{File} = basename($F{File});
      $readMetrics{$F{Sample}} = \%F;
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
  # If we don't have assembly metrics from the old version of assembleAll.pl, then
  # maybe we're looking for quast metrics.
  else {
    for my $samplename(keys(%$sampleInfo)){
      my $quastReport = "$dir/SneakerNet/assemblies/$samplename/quast/report.tsv";
      next if(!-e $quastReport);

      open(my $fh, $quastReport) or die "ERROR: could not read $quastReport: $!";
      while(<$fh>){
        chomp;
        my ($key, $value) = split(/\t/, $_);

        # I don't really need a pound sign for number sign.
        $key =~ s/#\s*//;

        # If the key describes a contig length >= than something,
        # then let's go with 1kb.
        if($key =~ /(.+?)\s*>=\s*(\d+)/){
          my $size = $2;
          next if($size != 1000);
          #$key = $1;
          $key =~ s/\s+\(.+//;
        }
        $key =~ s/\s+/_/g;

        $assemblyMetrics{$samplename}{$key} = $value;
      }
      close $fh;

      # change some headers to my own legacy names
      $assemblyMetrics{$samplename}{largestContig} = $assemblyMetrics{$samplename}{Largest_contig};
      $assemblyMetrics{$samplename}{genomeLength} = $assemblyMetrics{$samplename}{Total_length};
      $assemblyMetrics{$samplename}{GC} = $assemblyMetrics{$samplename}{'GC_(%)'};
      $assemblyMetrics{$samplename}{numContigs} = $assemblyMetrics{$samplename}{contigs};
      $assemblyMetrics{$samplename}{percentNs} = $assemblyMetrics{$samplename}{'N\'s_per_100_kbp'} / 10e5;
      # still don't have: avgContigLength assemblyScore effectiveCoverage expectedGenomeLength 
      #                   kmer21 minContigLength 
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
    if($samplename=~/^(Undetermined)/i){
      logmsg "Sample name $samplename contains '$1' and will therefore be skipped." if($$settings{debug});
      next;
    }

    # Get the name of all files linked to this file through the sample.
    $$sampleInfo{$samplename}{fastq}//=[]; # Set {fastq} to an empty list if it does not exist

    # Get metrics from the fastq files
    my $totalCoverage = 0;

    # Uncover the total coverage for the sample
    if(defined $readMetrics{$samplename}{coverage}){
      $totalCoverage = $readMetrics{$samplename}{coverage};
      logmsg "Coverage calculation: Sample $samplename => ${totalCoverage}x" if($$settings{debug});
      if($totalCoverage < $$sampleInfo{$samplename}{taxonRules}{coverage}){
        logmsg "  I will fail this sample for coverage: $samplename" if($$settings{debug});
        $fail{coverage} = 1;
      } else {
        logmsg "  I will not fail this sample for coverage: $samplename" if($$settings{debug});
        $fail{coverage} = 0;
      }
    }

    # Check R1 and R2 quality to be above threshold
    $$sampleInfo{$samplename}{taxonRules}{quality} ||= 0;
    # Quality scores are either defined or are the biggest int possible ~0.
    my $Q1 = $readMetrics{$samplename}{avgQuality1} // ~0;
    my $Q2 = $readMetrics{$samplename}{avgQuality2} // ~0;
    my $qThreshold = $$sampleInfo{$samplename}{taxonRules}{quality};
    if($Q1 >= $qThreshold && $Q2 >= $qThreshold){
      logmsg "Sample passed quality: $samplename" if($$settings{debug});
      $fail{quality} = 0;
    } else {
      logmsg "Sample failed quality: $samplename (Q1:$Q1 Q2:$Q2)" if($$settings{debug});
      $fail{quality} = 1;
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
        logmsg "Percent Ns: I cannot tell if I need to fail $samplename" if($$settings{debug});
        $failAssembly{minNs} = -1;
      }
      elsif($percentNs > $minNs){
        logmsg "Percent Ns: I will fail $samplename because it has $percentNs percent Ns" if($$settings{debug});
        $failAssembly{minNs} = 1;
      }
      elsif($percentNs <= $minNs){
        logmsg "Percent Ns: I will not fail $samplename because it has $percentNs percent Ns" if($$settings{debug});
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
  --debug    Explain why things fail
";
  exit(0);
}

