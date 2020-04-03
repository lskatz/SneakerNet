#!/usr/bin/env perl
# Use Kraken for shotgun metagenomics

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/cp mv/;
use File::Temp qw/tempdir/;
use File::Find qw/find/;
use File::Spec::Functions qw/abs2rel rel2abs/;
use FindBin;

use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readConfig samplesheetInfo_tsv command logmsg/;

our $VERSION = "1.0";
our $CITATION= "Kraken plugin for metagenomics by Lee Katz.";

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
      'ktImportKrona (Krona)'     => 'ktImportKrona | grep "/" | grep -P -m 1 -o "KronaTools .*ktImportKrona"',
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=tempdir("$0XXXXXX",TMPDIR=>1, CLEANUP=>1);

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet/forEmail";
  mkdir "$dir/SneakerNet/kraken";
  my $outHtml = "$dir/SneakerNet/forEmail/kraken.metagenomics.html";

  # Get all HTML files.
  # I didn't use glob here to stay safe, in case no
  # sample subdirectories are not there.  Or any other
  # unforseeable things.
  my @kronaHtml;
  find({no_chdir=>1,  wanted=>sub{
    return if(!/report\.html$/);
    my $sampleName = basename($File::Find::dir);
    my $target = "$$settings{tempdir}/$sampleName.html";
    logmsg "cp $_ => $target";
    cp($_, $target);

    push(@kronaHtml, $target);
  }}, ("$dir/SneakerNet/kraken"));
  my $kronaPosArgs = join(" ", @kronaHtml);
  logmsg $kronaPosArgs;
  
  unlink($outHtml);
  system("ktImportKrona $kronaPosArgs -o $outHtml");
  die "ERROR with ktImportKrona\n  ".join("\n  ",@kronaHtml) if $?;
  
  recordProperties($dir,{version=>$VERSION,krakenDatabase=>$$settings{KRAKEN_DEFAULT_DB},html=>$outHtml,table=>"$dir/SneakerNet/forEmail/kraken.metagenomics.tsv"});

  return 0;
}

sub usage{
  print "Summarize metagenomics analysis from kraken on a run
  Usage: $0 MiSeq_run_dir
  --numcpus 1
  --version
  --force         Overwrite any previous results
";
  exit(0);
}

