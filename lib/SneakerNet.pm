package SneakerNet;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Data::Dumper;

use FindBin qw/$Bin $Script $RealBin $RealScript/;

our @EXPORT_OK = qw(
  readConfig samplesheetInfo 
  command logmsg
);

sub logmsg{print STDERR "$0: @_\n";}

sub readConfig{
  #my @file=glob("$FindBin::RealBin/../config/*");
  my @file=glob(dirname($INC{'SneakerNet.pm'})."/../config/*");
  my $settings={};
  for(@file){
    open(CONFIGFILE,$_) or die "ERROR: could not open config file $_: $!";
    my $key=basename $_;
    while(<CONFIGFILE>){
      s/^\s+|\s+$//g; # trim
      next if(/^$/);
      next if(/^#/);
      my $configLine=[split(/\t/,$_)];
      push(@{ $$settings{$key} },$configLine);
    }
    close CONFIGFILE;
  }
  return $settings;
}

sub samplesheetInfo{
  my($samplesheet,$settings)=@_;

  my $config=readConfig();

  my $section="";
  my @header=();
  my %sample;
  open(SAMPLE,$samplesheet) or die "ERROR: could not open sample spreadsheet $samplesheet: $!";
  while(<SAMPLE>){
    s/^\s+|\s+$//g; # trim whitespace

    if(/^\[(\w+)\]/){  # [sectionname]
      $section=lc($1);
      my $header=<SAMPLE>;
      $header=~s/^\s+|\s+$//g; # trim whitespace
      @header=split(/,/,lc($header));
      next;
    }
    if($section eq "data"){
      my %F;
      @F{@header}=split(/,/,$_);
      for my $keyvalue(split(/;/,lc($F{description}))){
        my($key,$value)=split(/=/,$keyvalue);
        $key=~s/^\s+|\s+$//g;      #whitespace trim
        $value=~s/^\s+|\s+$//g;    #whitespace trim
        #$F{$key}={} if(!$F{$key});
        #$F{$key}{$value}++;
        if($F{$key}){
          if(ref($F{$key}) ne 'ARRAY'){
            $F{$key}=[$F{$key}];
          }
          push(@{ $F{$key} }, $value);
        } else {
          $F{$key}=$value;
        }
      }
      delete($F{description});

      # What taxon is this if not listed?
      if(!$F{species}){
        for my $taxonArr(@{ $$config{genomeSizes} }){
          my($regex,$size,$possibleTaxon)=@$taxonArr;
          if($F{sample_id}=~/$regex/){
            $F{species}=$possibleTaxon;
            last;
          }
        }
      }

      # Route needs to be an array.
      if(!defined($F{route})){
        $F{route}=[];
      } elsif(ref($F{route}) ne 'ARRAY'){
        $F{route}=[];
      }
      
      $sample{$F{sample_id}}=\%F;
    }
  }
  close SAMPLE;

  # Try to associate samples to files
  # Warning: this adds a mix of strings into a set of hashes and so
  # the variable type (ref) needs to be checked sometimes.
  my %fastqToName;
  while(my($samplename,$sampleinfo)=each(%sample)){
    my @possibleFastq=glob(dirname($samplesheet)."/$samplename*.fastq.gz");
    $sample{$samplename}{fastq}=\@possibleFastq;
    
    # Make some links from file to sample
    for my $fastq(@possibleFastq){
      $fastqToName{$fastq}=$samplename;
    }
  }
  %sample=(%sample,%fastqToName);

  return \%sample;
}

sub command{
  my($command,$settings)=@_;
  logmsg "COMMAND\n  $command" if($$settings{debug});
  system($command);
  if($?){
    my $msg="ERROR running command\n  $command";
  }
}


