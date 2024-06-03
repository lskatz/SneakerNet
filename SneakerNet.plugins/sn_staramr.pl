#!/usr/bin/env perl
# Hello World example

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv cp/;
use Cwd qw/realpath/;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "2.0";
our $CITATION = "StarAMR plugin by Lee Katz and Jess Chen";

my($basename, $thisDir) = fileparse $0;
local $0=$basename;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  my @exe = qw(staramr blastn mlst);
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      exe => \@exe,
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
  if(! -d "$dir/SneakerNet/staramr"){
    mkdir "$dir/SneakerNet/staramr";
  }

  staramr($dir, $settings);

  # Output file headers
  my @header = ("Sample", "Assembly", "Genotype", "Predicted Phenotype");
  my $outfile = "$dir/SneakerNet/forEmail/staramr.tsv";

  open(my $fh, '>', $outfile) or die "ERROR: could not write to $outfile: $!";
  print $fh join("\t", @header)."\n";
  for my $sampleDir(glob("$dir/SneakerNet/staramr/*")){
    next if(! -d $sampleDir);
    open(my $summaryFh, '<', "$sampleDir/summary.tsv") or die "ERROR: could not read from $sampleDir/summary.tsv: $!";
    my $header = <$summaryFh>;
    chomp($header);
    my @thisHeader = split(/\t/, $header);

    # Get the result
    my $res = <$summaryFh>;
    chomp($res);
    my @F = split(/\t/, $res);

    # Close out the file
    close $summaryFh;

    # Label the results
    my %F;
    @F{@thisHeader} = @F;
    # Extra labels
    $F{Sample} = basename($sampleDir);
    $F{Assembly} = $F{"Isolate ID"};

    # Print to a combined table
    for my $h(@header){
      print $fh "$F{$h}\t";
    }
    print $fh "\n";
  }
  print $fh "# The predicted phenotypes/drug resistances are for microbiological resistance and not clinical resistance. Predictions are provided with support from the NARMS/CIPARS Molecular Working Group with an emphasis on Salmonella, Shigella, E. coli, and Campylobacter and are continually being improved.\n";
  close $fh;

  my $rawMultiQC = makeMultiQC($dir, $settings);

  recordProperties($dir,{exe=>\@exe,version=>$VERSION, table=>$outfile, mqc=>$rawMultiQC});

  # staramr has additional properties for its database
  logmsg "Recording database properties";
  open(my $dbPropFh, " staramr db info | ") or die "ERROR: could not get properties for staramr database: $!";
  while(<$dbPropFh>){
    chomp;
    if(/(\S+)\s+(=?)\s*(.+)/){
      my ($key,$value) = ($1, $3);
      recordProperties($dir,{$key=>$value});
      logmsg "$key .. $value";
    }
  }
  close $dbPropFh;

  logmsg "Output table is in $outfile";

  return 0;
}

# Make a table suitable for MultiQC
# Example found at https://github.com/MultiQC/test-data/blob/main/data/custom_content/issue_1883/4056145068.variant_counts_mqc.tsv
sub makeMultiQC{
  my($dir, $settings) = @_;
  my $intable = "$dir/SneakerNet/forEmail/staramr.tsv";
  my $mqcDir  = "$dir/SneakerNet/MultiQC-build";
  my $outtable= "$mqcDir/staramr_mqc.tsv";
  mkdir($mqcDir);

  my $plugin = basename($0);
  my $anchor = basename($0, ".pl");

  my $docLink = "<a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>&#128196;</a>";
  my $pluginLink = "<a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'><span style='font-family:monospace;font-size:small'>1011</span></a>";

  open(my $outFh, ">", $outtable) or die "ERROR: could not write to multiqc table $outtable: $!";
  print $outFh "#id: $anchor'\n";
  print $outFh "#section_name: \"AMR\"\n";
  print $outFh "#description: \"Staramr AMR results<br />$plugin v$VERSION $docLink $pluginLink\"\n";
  print $outFh "#anchor: '$anchor'\n";
  # Print the rest of the table
  open(my $fh, $intable) or die "ERROR: could not read table $intable: $!";
  while(<$fh>){
    next if(/^#/);
    s/ /_/g;
    print $outFh $_;
  }
  close $fh;

  return $outtable;
}

sub staramr{
  my($dir, $settings)=@_;

  my $samples = samplesheetInfo_tsv("$dir/samples.tsv", $settings);

  while(my($sampleName, $s) = each(%$samples)){
    # Jess's shell script:
    # detectARDs.sh genus assembly R1 R2 strainID
    # staramr search --pointfinder-organism $1 --exclude-genes-file $1_genes_to_exclude.tsv --pid-threshold 90 --percent-length-overlap-resfinder 50  -o $5.staramr $2

    # Specify the output directory in temporary land
    my $tempdir = "$$settings{tempdir}/$sampleName";
    system("rm -rf $tempdir"); # ensure this directory does not exist yet
    # ... and the final directory
    my $outdir = "$dir/SneakerNet/staramr/$sampleName";

    if(!$$settings{force} && -d $outdir){
      logmsg "Already found $outdir. Skipping";
      next;
    }
    system("rm -rf $outdir");

    # Get the genome assembly
    my @asm = glob("$dir/SneakerNet/assemblies/$sampleName/*.fasta");
    my $asm = $asm[0];

    if(!defined($asm) || !-e $asm || -s $asm < 30000){
      logmsg "Assembly for $sampleName is too small. Skipping.";
      next;
    }

    # Run staramr
    logmsg "staramr on $sampleName";
    my $staramrXopts = "";
    if(my $pointfinderOrganism = $$s{taxonRules}{pointfinder}){
      logmsg "Pointfinder organism found for $sampleName: $pointfinderOrganism";
      $staramrXopts .= "--pointfinder-organism $pointfinderOrganism ";
      
      # But what if the installation directory doesn't have the right
      # database? See if we have the database locally in
      # Sneakernet/db/staramr
      my @possibleDir = ("$thisDir/../db/staramr/pointfinder/$pointfinderOrganism",
                         "$thisDir/../db/staramr/update/pointfinder/$pointfinderOrganism",
                        );

      for my $possibleDir(@possibleDir){
        if(-e $possibleDir){
          $possibleDir =~ s/staramr.*/staramr/;
          my $path = realpath("$possibleDir");
          logmsg "Using pointfinder database $possibleDir";
          $staramrXopts .= "--database $path";
          last;
        } 
        #else {logmsg "NOT FOUND $possibleDir";}
      }
    } else {
      logmsg "Pointfinder organism not found for $sampleName";
    }
    my $staramrLog = "$$settings{tempdir}/$sampleName.staramr.log";
    my $command = "staramr search $staramrXopts --pid-threshold 90 --percent-length-overlap-resfinder 50 --output-dir $tempdir $asm 2>$staramrLog";
    system($command);
    if($?){
      die "ERROR running\n  $command\nERROR was\n".`cat $staramrLog`;
    }

    system("mv $tempdir $outdir");
  }

  return 1;
}
    

sub usage{
  print "Run StarAMR resistance finding

  Usage: $0 MiSeq_run_dir
  ";
  exit 0;
}

