#!/usr/bin/env perl
# Simply run kraken so that results are available to other plugins.
# Do not make any high-level conclusions.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/cp mv/;
use File::Temp qw/tempdir/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readTsv exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "1.1";
our $CITATION= "Kraken plugin by Lee Katz.  Uses Kraken1.";

my %errors = ();

# Get the executable directories
my $tmpSettings=readConfig();

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help debug tempdir=s numcpus=i force)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
      zip       => 'zip --version | grep "This is Zip"',
      kraken    => 'kraken --version | grep -m 1 version',
      'kraken-translate (Kraken)' => 'kraken-translate --version | grep -m 1 version',
      'kraken-report (Kraken)'    => 'kraken-report --version | grep -m 1 version',
      'ktImportText (Krona)'     => 'ktImportText | grep "/" | grep -P -m 1 -o "KronaTools .*ktImportText"',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{KRAKEN_DEFAULT_DB} ||= die "ERROR: KRAKEN_DEFAULT_DB needs to be defined under config/settings.conf";
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/forEmail";
  mkdir "$dir/SneakerNet/kraken";

  # Pass read metrics along to the rest of the analysis through $$settings{readMetrics}
  my $readMetrics = readTsv("$dir/SneakerNet/forEmail/readMetrics.tsv");
  while(my($fastq, $metrics) = each(%$readMetrics)){
    # Save the full path
    my $fastqPath = File::Spec->rel2abs($fastq,$dir);
    $$settings{readMetrics}{$fastqPath} = $metrics;
  }
  
  my $outdir=runKrakenOnDir($dir,$settings);

  my $errorsMsg = join(" ", keys(%errors));
  recordProperties($dir,{version=>$VERSION,krakenDatabase=>$$settings{KRAKEN_DEFAULT_DB}, errors=>$errorsMsg,});

  return 0;
}

sub runKrakenOnDir{
  my($dir,$settings)=@_;

  # Expect output to go here
  my $outdir="$dir/SneakerNet/kraken";
  # Sanity check
  if(!-e $outdir){
    die "ERROR: output directory does not exist $outdir";
  }

  my $sampleInfo=samplesheetInfo_tsv("$dir/samples.tsv",$settings);

  my @report; # reporting contamination in an array, in case I want to sort it later
  while(my($sampleName,$s)=each(%$sampleInfo)){

    my $sampledir="$outdir/$sampleName";
    mkdir $sampledir;

    logmsg "Running Kraken on $sampleName";
    logmsg "  Database: $$settings{KRAKEN_DEFAULT_DB}";
    my $krakenWorked=runKraken($s,$sampledir,$settings);

    # Check for anything odd in the report
    if(-e "$sampledir/kraken.report"){
      open(my $fh, "<", "$sampledir/kraken.report") or logmsg "Could not open $sampledir/kraken.report: $!";
      while(<$fh>){
        # Just look at the first integer and ignore decimals
        # for the sake of simplicity
        if(/(\d+)/){
          my $percentMatch = $1;
          if($percentMatch > 100){
            $errors{"Internal error: There is at least one sample that has at least 100% of the hits"}++;
            last;
          }
        }
      }
      close $fh;
    }

    if(!$krakenWorked){
      logmsg "Kraken was not completed successfully on $sampleName. I will not output results for $sampleName";
      $errors{"Kraken did not complete successfully for at least one sample"}++;
      next;
    }

  }

  return $outdir;
}

sub runKraken{
  my($sample,$sampledir,$settings)=@_;
  my $sampleName = $$sample{sample_id};

  my $html="$sampledir/report.html";
  if(-e $html && !$$settings{force}){
    return 1;
  }
  
  # Skip any samples without reads or assemblies, ie, samples that are misnamed or not sequenced.
  # There is no way to predict how a sample is misnamed and so it does not fall under
  # this script's purview.
  my $inputType = "READS";
  if(!defined($$sample{fastq}) || !@{ $$sample{fastq} }){
    logmsg "WARNING: I could not find the reads for $sampleName .";
    $inputType = "";
    return 0;
    if(!defined($$sample{asm}) || (ref($$sample{asm} eq 'ARRAY') &&  !@{ $$sample{asm} })){
      logmsg "WARNING: I could not find the assembly for $sampleName .";
      return 0; # no reads or asm -- give up and return 0
    } else {
      $inputType = "ASM";
    }
  }

  # Force an array
  if(ref($$sample{fastq}) ne 'ARRAY'){
    $$sample{fastq} = [$$sample{fastq}];
  }
  my @fastq = @{$$sample{fastq}};

  my $asm = $$sample{asm};

  if(@fastq < 2 && -e $asm){
    return runKrakenAsm(@_);
  }
  if(@fastq == 1){
    return runKrakenSE(@_);
  }
  elsif(@fastq == 2){
    return runKrakenPE(@_);
  } 
  elsif(@fastq > 2) {
    my %sampleCopy = %$sample;
    splice(@{ $sampleCopy{fastq} }, 2);
    return runKrakenPE(\%sampleCopy, $sampledir, $settings);
  }

  die "INTERNAL ERROR";
  return 0;
}
sub runKrakenAsm{
  my($sample,$sampledir,$settings)=@_;

  my $sampleName = $$sample{sample_id};
  my $asm = $$sample{asm};

  # Skip small file sizes.
  # TODO: use something better like readMetrics.pl 
  if(-s $asm < 1000){
    logmsg "The assembly is too small for $sampleName. Skipping";
    return 0;
  }

  my $html="$sampledir/report.html";

  # Run basic kraken command
  command("kraken --fasta-input $asm --db=$$settings{KRAKEN_DEFAULT_DB} --threads $$settings{numcpus} --output $sampledir/kraken.out.tmp ");
  mv("$sampledir/kraken.out.tmp", "$sampledir/kraken.out");

  # Create the taxonomy but normalize for contig length
  # I stole my own code from https://github.com/lskatz/lskScripts/blob/master/scripts/translate-kraken-contigs.pl
  my %length;     # Contig lengths
  my %percentage; # Percentage of all nucleotides
  open(my $fh, '<', "$sampledir/kraken.out") or die "ERROR: could not read $sampledir/kraken.out: $!";
  while(<$fh>){
    chomp;
    my($classified,$seqname,$taxid,$length,$kmerTaxid)=split(/\t/,$_);
    if($classified eq 'U'){
      $percentage{'unclassified'}+=$length;
    } else {
      $length{$seqname} = $length;
    }
  }
  close $fh;

  # kraken-translate but tally all the sequence lengths
  open(my $translateFh, "kraken-translate $sampledir/kraken.out | ") or die "ERROR: could not run kraken-translate on $sampledir/kraken.out: $!";
  while(<$translateFh>){
    chomp;
    my($seqname,$taxonomyString)=split(/\t/,$_);
    $taxonomyString=~s/\s+/_/g;
    $taxonomyString=~s/;/\t/g;
    $percentage{$taxonomyString}+=$length{$seqname};
  }
  close $translateFh;

  # Create the kraken-translate file
  my $taxonomyFile = "$sampledir/kraken.taxonomy";
  open(my $taxonomyFh, '>', "$taxonomyFile.tmp") or die "ERROR: could not write to $taxonomyFile.tmp: $!";
  while(my($taxonomyString,$sliceOfPie)=each(%percentage)){
    print $taxonomyFh join("\t",$sliceOfPie,$taxonomyString)."\n";
  }
  close $taxonomyFh;
  mv("$taxonomyFile.tmp", $taxonomyFile) or die $!;

  command("kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");

  filterReport($sampledir, $settings);

  command("ktImportText -o $html $sampledir/kraken.taxonomy,$sampleName");

  unlink("$sampledir/kraken.out");
  
  return 1;
}


sub runKrakenSE{
  my($sample,$sampledir,$settings)=@_;

  my $sampleName = $$sample{sample_id};
  
  # Skip small file sizes.
  # TODO: use something better like readMetrics.pl 
  for(@{ $$sample{fastq} }){
    if(-s $_ < 10000){
      logmsg "There are few reads in $sampleName. Skipping.";
      return 0;
    }
  }
  
  my $html="$sampledir/report.html";

  my $reads = $$sample{fastq}[0];

  return 0 if(!$reads);

  command("kraken --fastq-input $reads --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out ");

  command("kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  command("kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");
  filterReport($sampledir, $settings);


  # To capture unclassified reads, we can get the third
  # column of the first row of the report file. This
  # information can be appended to the taxonomy file
  # on the last line.
  open(my $reportFh, "<", "$sampledir/kraken.report") or die "ERROR: could not read $sampledir/kraken.report: $!";
  my $firstLine=<$reportFh>;
  close $reportFh;
  my $unclassifiedReadsCount=(split(/\t/, $firstLine))[2];
  open(my $taxFh, ">>", "$sampledir/kraken.taxonomy") or die "ERROR: could not append to $sampledir/kraken.taxonomy: $!";
  print $taxFh $unclassifiedReadsCount."\n";
  close $taxFh;

  command("ktImportText -o $html $sampledir/kraken.taxonomy,$sampleName");

  # Go ahead and remove kraken.out which is a huge file
  unlink("$sampledir/kraken.out");

  if(! -e "$sampledir/kraken.taxonomy"){
    return 0;
  }

  return 1;
}

sub runKrakenPE{
  my($sample,$sampledir,$settings)=@_;

  my $sampleName = $$sample{sample_id};
  
  # Skip files with few reads.
  for(@{ $$sample{fastq} }){
    my $fastqFilename = basename($_);
    my $numReads = $$settings{readMetrics}{$fastqFilename}{numReads};
    if(defined($numReads) && $numReads < 100){
      logmsg "There are few reads in $sampleName (N = $numReads). Skipping.";
      return 0;
    }
  }
  
  my $html="$sampledir/report.html";

  my @twoReads = (@{$$sample{fastq}})[0,1];
  my $reads="'".join("' '", @twoReads)."'";
  return 0 if(!$reads);
  
  command("kraken --fastq-input --paired $reads --db=$$settings{KRAKEN_DEFAULT_DB} --gzip-compressed --quick --threads $$settings{numcpus} --output $sampledir/kraken.out ");

  command("kraken-translate --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out | cut -f 2- | sort | uniq -c | perl -lane '
    s/^ +//;   # remove leading spaces
    s/ +/\t/;  # change first set of spaces from uniq -c to a tab
    s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
    print;
    ' | sort -k1,1nr > $sampledir/kraken.taxonomy
  ");

  command("kraken-report --db $$settings{KRAKEN_DEFAULT_DB} $sampledir/kraken.out > $sampledir/kraken.report");
  filterReport($sampledir, $settings);


  # To capture unclassified reads, we can get the third
  # column of the first row of the report file. This
  # information can be appended to the taxonomy file
  # on the last line.
  open(my $reportFh, "<", "$sampledir/kraken.report") or die "ERROR: could not read $sampledir/kraken.report: $!";
  my $firstLine=<$reportFh>;
  close $reportFh;
  my $unclassifiedReadsCount=(split(/\t/, $firstLine))[2];
  open(my $taxFh, ">>", "$sampledir/kraken.taxonomy") or die "ERROR: could not append to $sampledir/kraken.taxonomy: $!";
  print $taxFh $unclassifiedReadsCount."\n";
  close $taxFh;

  command("ktImportText -o $html $sampledir/kraken.taxonomy,$sampleName");

  # Go ahead and remove kraken.out which is a huge file
  unlink("$sampledir/kraken.out");

  if(! -e "$sampledir/kraken.taxonomy"){
    return 0;
  }

  return 1;
}

sub filterReport{
  my($sampledir, $settings) = @_;

  # Also filter the kraken report to remove small percentages
  # in pure perl to make it more portable
  open(my $reportFh, "<", "$sampledir/kraken.report") or die "ERROR: could not open $sampledir/kraken.report: $!";
  open(my $filteredFh, ">", "$sampledir/kraken.filtered.report") or die "ERROR: could not write to $sampledir/kraken.filtered.report: $!";
  while(<$reportFh>){
    # Just check out the first digits and don't sweat the decimals.
    if(/(\d+)/){
      my $percentageInt = $1;
      next if($percentageInt < 1);

      print $filteredFh $_;
    }
  }
  close $reportFh;
  close $filteredFh;

  return "$sampledir/kraken.filtered.report";
}

sub usage{
  print "Runs kraken in SneakerNet. Results are made available to other plugins.
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
  --force         Overwrite any previous results
";
  exit(0);
}

