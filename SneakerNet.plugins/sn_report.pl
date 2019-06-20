#!/usr/bin/env perl
# Creates a report for SneakerNet

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use File::Temp qw/tempdir/;
use File::Spec;
use Cwd qw/realpath/;
use POSIX qw/strftime/;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/recordProperties readProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec/;

our $VERSION = "1.3";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version help force tempdir=s debug numcpus=i)) or die $!;
  if($$settings{version}){
    print $VERSION."\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/forEmail";

  # Special to this plugin: record any properties before
  # reading the properties.
  recordProperties($dir,{
    version=>$VERSION,
    date=>strftime("%Y-%m-%d", localtime()),
    time=>strftime("%H:%M:%S", localtime()),
  });

  my $properties = readProperties($dir);

  my $html="";

  $html .= htmlHeaders();
  $html .= "<H1>QC report for ".basename(realpath($dir))."</H1>\n";
  $html .= "<p class='genericInfo' style='margin-bottom:1em'>SneakerNet reporting plugin version $VERSION</p>\n";
  $html .= "<H2>Table of Contents</H2>\n";
  $html .= "<ul>\n";
  for my $plugin(sort keys(%$properties)){
    $html .= "<li><a href='#plugin-$plugin'>$plugin</a> v$$properties{$plugin}{version}</li>\n"; # TODO add version here too
  }
  $html .= "</ul>\n";

  for my $plugin(sort keys(%$properties)){
    $html .= "<div class='pluginSplash'>\n";
    $html .= "<H2 id='plugin-$plugin'>$plugin</H2>\n";
    $html .= report($plugin, $properties, $settings);
    $html .= "</div>\n";
  }
  $html .= htmlFooters();

  my $outfile = "$dir/SneakerNet/forEmail/report.html";
  open(my $fh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  print $fh $html;
  close $fh;

  return 0;
}

sub report{
  my($plugin, $p, $settings) = @_;

  my $html = "";

  my @sortedKeys = sort{
    if($a =~ /version/){
      return 1;
    }
    if($b =~ /version/){
      return -1;
    }

    return $a cmp $b;
  } keys(%{ $$p{$plugin} });

  for my $key(@sortedKeys){
    my $value = $$p{$plugin}{$key};
    if($value =~ /\.(tsv|csv)$/i){
      my $type = $1;
      $html .= tableHtml($plugin, $key, $value, $type);
    }
    else{
      $html .= genericHtml($plugin, $key, $value);
    }
  }

  return $html;
}

sub tableHtml{
  my($plugin, $key, $value, $type) = @_;

  # Set the separator
  my $sep = "\t";
  if($type eq 'csv'){
    $sep = ',';
  }

  my $html = "";

  $html.="<table cellspacing=0>\n";

  $html.="<tbody>\n";
  
  my $numColumns= 1;
  my $seen_footer = 0;
  my $rowCounter = 0;
  open(my $fh, $value) or die "ERROR: could not open $value: $!";
  while(<$fh>){
    chomp;

    if(/^#|^NOTE/i){
      $seen_footer = 1;
      s/^#//;
      $html .= "</tbody>\n<tfoot>\n";
    }
    if($rowCounter == 0){
      $html.="<thead>\n";
    }

    my @F = split /$sep/;
    my $colspan = 2;
    if($seen_footer){
      $colspan = $numColumns;
    } else {
      $numColumns = scalar(@F);
      $colspan = 1;
    }

    $html .= "<tr>\n";
    for my $value(@F){
      $html .= "<td colspan='$colspan' seen_footer='$seen_footer'>$value</td>"
    }
    $html .= "\n</tr>\n";
    if($seen_footer){
      $html .= "<tfoot>\n";
    }
    if($rowCounter == 0){
      $html.="</thead>\n";
    }

    $rowCounter++;
  }
  close $fh;

  $html.="</table>\n";

  return $html;
}
sub genericHtml{
  my($plugin, $key, $value) = @_;
  my $html = "";

  my $class="genericInfo";
  if($key =~ /version/i){
    $class.=" version";
  }
  $html .= "<p class='$class'>$key: $value</p>\n";
  return $html;
}

sub htmlHeaders{
  my($settings)=@_;
  my $html = "";

  $html .= "<html><head><title>SneakerNet report</title>\n";

  $html .= "<style>\n";
  $html .= "h1    {font-weight:bold; font-size:24px; margin:3px 0px;}\n";
  $html .= "h2    {font-weight:bold; font-size:16px; margin:3px 0px;}\n";
  $html .= "body  {font-size:12px; font-family:-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif, Apple Color Emoji, Segoe UI Emoji, Segoe UI Symbol;}\n";
  $html .= "table {border:1px solid black;}\n";
  $html .= "td    {border:1px solid #999999; margin:0px;}\n";
  $html .= "th    {border:1px solid #009900; margin:0px;}\n";
  $html .= "thead {font-weight:bold;}\n";
  $html .= "tbody {color:black;} \n";
  $html .= "tfoot {font-size:75%; color:#DD3333;}\n";
  $html .= ".genericInfo {background-color:#EEEEEE; border: 1px solid #666666; margin:2px 0px; padding:1px;}\n";
  $html .= ".version {font-size:75%; padding:1px;font-family:monospace; font-size:10px;}\n";
  $html .= ".pluginSplash{background-color:#FAFFFF;border:1px solid black; margin:6px 0px; padding:1px;}\n";
  $html .= "</style>\n";

  $html .= "</head>\n";
  
  return $html;
}

sub htmlFooters{
  my($settings)=@_;

  my $html = "";

  $html .= "</html>\n";

  return $html;
}

sub usage{
  "Make an HTML report from SneakerNet results
  Usage: $0 MiSeq_run_dir
  --version
  "
}

