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

our $VERSION = "1.4";

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
  $html .= "<p class='genericInfo'>SneakerNet version $SneakerNet::VERSION</p>\n";
  #$html .= "<H2>Table of Contents</H2>\n";
  #$html .= "<ul>\n";
  #for my $plugin(sort keys(%$properties)){
  #  $html .= "<li><a href='#plugin-$plugin'>$plugin</a> v$$properties{$plugin}{version}</li>\n"; # TODO add version here too
  #}
  #$html .= "</ul>\n";

  for my $plugin(sort keys(%$properties)){
    my $inputID = "menu-$plugin";
       $inputID =~s/\.+/-/g;
    $html .= "<div class='pluginSplash'>\n";
    $html .= "<input type='checkbox' class='collapsible' id='$inputID' />\n";
    $html .= "<label class='collapsible' for='$inputID'>$plugin v$$properties{$plugin}{version}</label>\n";
    $html .= "<div class='pluginContent'>\n";
    #$html .= "<H2 id='plugin-$plugin'>$plugin</H2>\n";
    $html .= report($plugin, $properties, $settings);
    $html .= "</div>\n";
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
  $html.='<meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">';
  #$html .= "<!DOCTYPE PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN' >\n";

  # Images encoded from https://www.iconfinder.com
  my $menuBase64 = 'data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiA/PjwhRE9DVFlQRSBzdmcgIFBVQkxJQyAnLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4nICAnaHR0cDovL3d3dy53My5vcmcvR3JhcGhpY3MvU1ZHLzEuMS9EVEQvc3ZnMTEuZHRkJz48c3ZnIGhlaWdodD0iMzJweCIgaWQ9IkxheWVyXzEiIHN0eWxlPSJlbmFibGUtYmFja2dyb3VuZDpuZXcgMCAwIDMyIDMyOyIgdmVyc2lvbj0iMS4xIiB2aWV3Qm94PSIwIDAgMzIgMzIiIHdpZHRoPSIzMnB4IiB4bWw6c3BhY2U9InByZXNlcnZlIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIj48cGF0aCBkPSJNNCwxMGgyNGMxLjEwNCwwLDItMC44OTYsMi0ycy0wLjg5Ni0yLTItMkg0QzIuODk2LDYsMiw2Ljg5NiwyLDhTMi44OTYsMTAsNCwxMHogTTI4LDE0SDRjLTEuMTA0LDAtMiwwLjg5Ni0yLDIgIHMwLjg5NiwyLDIsMmgyNGMxLjEwNCwwLDItMC44OTYsMi0yUzI5LjEwNCwxNCwyOCwxNHogTTI4LDIySDRjLTEuMTA0LDAtMiwwLjg5Ni0yLDJzMC44OTYsMiwyLDJoMjRjMS4xMDQsMCwyLTAuODk2LDItMiAgUzI5LjEwNCwyMiwyOCwyMnoiLz48L3N2Zz4=';
  my $closedMenuBase64 = 'data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiA/PjxzdmcgZmlsbD0ibm9uZSIgaGVpZ2h0PSIyNCIgc3Ryb2tlPSIjMDAwIiBzdHJva2UtbGluZWNhcD0icm91bmQiIHN0cm9rZS1saW5lam9pbj0icm91bmQiIHN0cm9rZS13aWR0aD0iMiIgdmlld0JveD0iMCAwIDI0IDI0IiB3aWR0aD0iMjQiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+PGxpbmUgeDE9IjE4IiB4Mj0iNiIgeTE9IjYiIHkyPSIxOCIvPjxsaW5lIHgxPSI2IiB4Mj0iMTgiIHkxPSI2IiB5Mj0iMTgiLz48L3N2Zz4=';

  $html .= "<html><head><title>SneakerNet report</title>\n";

  $html .= "<style>\n";
  $html .= "h1    {font-weight:bold; font-size:24px; margin:3px 0px;}\n";
  $html .= "h2    {font-weight:bold; font-size:16px; margin:3px 0px;}\n";
  $html .= "body  {font-size:12px; font-family:-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif, Apple Color Emoji, Segoe UI Emoji, Segoe UI Symbol;}\n";
  $html .= "table {border:1px solid black;}\n";
  $html .= "td    {border:1px solid #999999; margin:0px; word-wrap: break-all;}\n";
  $html .= "th    {border:1px solid #009900; margin:0px; word-wrap: break-all;}\n";
  $html .= "thead {font-weight:bold; background-color:#BBE;}\n";
  $html .= "tbody {color:black;} \n";
  $html .= "tfoot {font-size:75%; color:#DD3333;}\n";
  $html .= ".genericInfo {background-color:#EEEEEE; border: 1px solid #666666; margin:2px 0px; padding:1px;}\n";
  $html .= ".genericInfo p {margin-bottom:1em;}\n";
  $html .= ".version {font-size:75%; padding:1px; margin-top:10px; font-family:monospace; font-size:10px;}\n";
  $html .= ".pluginSplash{background-color:#FAFFFF;border:1px solid black; margin:6px 0px; padding:1px;}\n";
  # https://codeburst.io/how-to-make-a-collapsible-menu-using-only-css-a1cd805b1390
  $html .= ".pluginContent{max-height:0px; overflow:hidden; }\n";
  $html .= "label.collapsible {
    display: block; 
    font-size:16px;
    color:black;
    font-weight:bold;
    cursor: pointer; 
    padding: 10px 0 10px 50px; 
    background-image: url(\"$menuBase64\");
    background-repeat: no-repeat;
    background-position-y: center;
    background-position-x: left;
  }\n";
  $html .= "input.collapsible {display:none;}\n";
  $html .= "input:checked ~ .pluginContent{max-height:100%; transition: all ease 1s;}\n";
  $html .= "input:checked ~ label {
    background-image: url(\"$closedMenuBase64\");
    transition: all ease 0.1s;
  }\n";
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

