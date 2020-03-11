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
use List::Util qw/min max/;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/exitOnSomeSneakernetOptions recordProperties readProperties readConfig samplesheetInfo_tsv command logmsg fullPathToExec passfail/;

our $VERSION = "2.2";
our $CITATION= "SneakerNet report by Lee Katz";

local $0=fileparse $0;
exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(version citation check-dependencies help force tempdir=s debug numcpus=i)) or die $!;
  exitOnSomeSneakernetOptions({
      _CITATION => $CITATION,
      _VERSION  => $VERSION,
    }, $settings,
  );

  usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{tempdir}||=File::Temp::tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);

  my $dir=$ARGV[0];
  mkdir "$dir/SneakerNet";
  mkdir "$dir/SneakerNet/forEmail";

  # Make the summary table before writing properties
  # so that it can get encoded with the rest of the 
  # tables in the html.
  my $file = "$dir/SneakerNet/forEmail/QC_summary.tsv";
  makeSummaryTable($file, $dir, $settings);

  my $pathFromDir = File::Spec->abs2rel(realpath($file), realpath($dir));

  # Special to this plugin: record any properties before
  # reading the properties.
  recordProperties($dir,{
    table => $pathFromDir,
    version=>$VERSION,
    date=>strftime("%Y-%m-%d", localtime()),
    time=>strftime("%H:%M:%S", localtime()),
    fullpath=>realpath($dir).'/SneakerNet',
  });

  my $properties = readProperties($dir);

  my $html="";

  $html .= htmlHeaders();
  $html .= "<H1>QC report for ".basename(realpath($dir))."</H1>\n";
  $html .= "<p class='genericInfo'><a href='https://github.com/lskatz/SneakerNet'>SneakerNet</a> version $SneakerNet::VERSION</p>\n";

  my @sortedPluginName = sort{
                        if($a eq basename($0)){
                          return -1;
                        }
                        elsif($b eq basename($0)){
                          return 1;
                        }
                        $a cmp $b;
                        } keys(%$properties);

  for my $plugin(@sortedPluginName){
    my $inputID = "menu-$plugin";
       $inputID =~s/\.+/-/g;

    # HTML attribute: should the checkbox be checked?
    # If so, the splash div will be open.
    my $checked = "";
    if($plugin eq basename($0)){
      $checked = 'checked';
    }
    my $version = $$properties{$plugin}{version};
    if(!$version){
      logmsg "WARNING: no version found for plugin $plugin. Setting to version '-1'";
      $version = -1;
    }

    $html .= "<div class='pluginSplash'>\n";
    $html .= "<input type='checkbox' class='collapsible' id='$inputID' $checked></input>\n";
    $html .= "<label class='collapsible' for='$inputID'>";
    $html .= "  $plugin v$version";
    $html .= "</label>\n";
    $html .= "<div class='pluginContent'>\n";
    $html .= report($dir, $plugin, $properties, $settings);
    $html .= "  <a title='documentation' href='https://github.com/lskatz/sneakernet/blob/master/docs/plugins/$plugin.md'>\n";
    $html .= "     <span style='font-family:monospace; border:solid pink 1px;margin:3px;'>&#128196;</span>\n"; # UTF-8 for page facing up
    #$html .= "    <img style='width:16px;height:16px;' alt='documentation' src='$documentationBase64' />\n";
    $html .= "  </a>";
    $html .= "  <a title='$plugin on github' href='https://github.com/lskatz/sneakernet/blob/master/SneakerNet.plugins/$plugin'>\n";
    $html .= "     <span style='font-family:monospace; border:solid pink 1px;margin:3px;'>1011</span>\n";
    #$html .= "    <img style='width:16px;height:16px;' alt='github' src='$githubBase64' />\n";
    $html .= "  </a>";
    $html .= "</div>\n"; # end div pluginContent
    $html .= "</div>\n"; # end div pluginSplash
  }
  $html .= htmlFooters($dir,$settings);

  my $outfile = "$dir/SneakerNet/forEmail/report.html";
  open(my $fh, ">", $outfile) or die "ERROR: could not write to $outfile: $!";
  print $fh $html;
  close $fh;
  logmsg "Report can be found in $outfile";

  return 0;
}

sub makeSummaryTable{
  my($outfile, $dir, $settings) = @_;

  # Gather some information
  my $sample   = samplesheetInfo_tsv("$dir/samples.tsv", $settings);
  my $passfail = passfail($dir, $settings);
  # Also add in contamination detection
  if(-e "$dir/SneakerNet/forEmail/kraken.tsv"){
    open(my $krakenFh, '<', "$dir/SneakerNet/forEmail/kraken.tsv") or logmsg "WARNING: kraken results were not found in $dir/SneakerNet/forEmail/kraken.tsv: $!";
    my $header = <$krakenFh>;
    chomp($header);
    my @header = split(/\t/, $header);
    while(<$krakenFh>){
      chomp;
      my @F = split(/\t/, $_);
      my %F;
      @F{@header} = @F;
      $F{PERCENTAGE_CONTAMINANT} //= 0;

      if($F{PERCENTAGE_CONTAMINANT} > 10){
        $$passfail{$F{NAME}}{kraken}=1;
      }
    }
    close $krakenFh;
  }
  
  # Read the readMetrics file into %readMetrics
  my %readMetrics = ();
  if(-e "$dir/readMetrics.tsv"){
    open(my $readMetricsFh,"$dir/readMetrics.tsv") or die "ERROR: could not open $dir/readMetrics.tsv: $!";
    my @rmHeader=split(/\t/,<$readMetricsFh>); chomp(@rmHeader);
    while(<$readMetricsFh>){
      chomp;
      my %F;
      @F{@rmHeader}=split(/\t/,$_);
      $F{File} = basename($F{File});
      $readMetrics{$F{File}} = \%F;
    }
    close $readMetricsFh;
  }

  my $happiness = $$settings{happiness_range} || [];
  if(ref($happiness) ne 'ARRAY' || !@$happiness){
    my @default = ('&#128515;', '&#128556;', '&#128561;');
    logmsg "WARNING: happiness_range was not set in settings.conf. Creating a default of @default";
    $happiness = \@default;
  }
  # Any knock on the perfection of a genome gets this penalty
  my $penalty = 1/scalar(@$happiness) * 100;

  my @tableRow; # to be sorted
  while(my($sampleName, $s) = each(%$sample)){
    my $score = 100;
    my $emojiIdx= 0;
    my @failure_code;
    my $failures = $$passfail{$sampleName};
    die "ERROR: sn_passfail.pl was not run on sample $sampleName" if(!$failures);
    while(my($failure_code, $is_failure) = each(%$failures)){
      if($is_failure==1){
        $score = $score - $penalty;
        $emojiIdx++;
        push(@failure_code, $failure_code);
      }
    }

    $score = sprintf("%0.0f", $score); # round
    $score = max($score, 0); # can't go negative on score
    
    # The emoji index is at most the last element of the emojis
    $emojiIdx = min($emojiIdx, scalar(@$happiness)-1);
    my $emoji = $$happiness[$emojiIdx];

    # Get the coverage and quality
    my $cov = "";
    my $qual= "";
    for my $fastq(@{ $$s{fastq} }){
      my $f = basename($fastq);
      $readMetrics{$f}{coverage}   //= -1;
      $readMetrics{$f}{avgQuality} //= -1;
      $cov .= sprintf("%0.0f",$readMetrics{$f}{coverage}).' ';
      $qual.= sprintf("%0.0f",$readMetrics{$f}{avgQuality}).' ';
    }
    
    @failure_code = ("None") if(!@failure_code);
    push(@tableRow, [$sampleName, $emoji, $score, join(", ", @failure_code),$qual, $cov, $$s{taxon}]);

  }
  my @sortedRow = sort{$$a[2] <=> $$b[2] || $$a[0] cmp $$b[0]} @tableRow;

  # Write the summary table
  open(my $fh, '>', $outfile) or die "ERROR: could not write to $outfile: $!";
  print $fh join("\t", qw(sample emoji score failure_code qual cov taxon))."\n";
  for my $row(@sortedRow){
    print $fh join("\t", @$row)."\n";
  }
  print $fh "# Scores start at 100 percent and receive a 33 percent penalty for each: low coverage, low quality, or high contamination in the Kraken report.\n";
  close $fh;
}

sub report{
  my($dir, $plugin, $p, $settings) = @_;

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
      $html .= tableHtml($dir, $plugin, $key, $value, $type);
    }
    else{
      $html .= genericHtml($plugin, $key, $value);
    }
  }

  return $html;
}

sub tableHtml{
  my($dir, $plugin, $key, $value, $type) = @_;

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
  my $abspath = File::Spec->rel2abs($value, $dir);
  open(my $fh, $abspath) or return "";
  #logmsg("WARNING: could not open $abspath: $!");
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
  #$html .= "<!DOCTYPE PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN' >\n";

  # Images encoded from https://www.iconfinder.com
  my $menuBase64 = 'data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiA/PjwhRE9DVFlQRSBzdmcgIFBVQkxJQyAnLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4nICAnaHR0cDovL3d3dy53My5vcmcvR3JhcGhpY3MvU1ZHLzEuMS9EVEQvc3ZnMTEuZHRkJz48c3ZnIGhlaWdodD0iMzJweCIgaWQ9IkxheWVyXzEiIHN0eWxlPSJlbmFibGUtYmFja2dyb3VuZDpuZXcgMCAwIDMyIDMyOyIgdmVyc2lvbj0iMS4xIiB2aWV3Qm94PSIwIDAgMzIgMzIiIHdpZHRoPSIzMnB4IiB4bWw6c3BhY2U9InByZXNlcnZlIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIj48cGF0aCBkPSJNNCwxMGgyNGMxLjEwNCwwLDItMC44OTYsMi0ycy0wLjg5Ni0yLTItMkg0QzIuODk2LDYsMiw2Ljg5NiwyLDhTMi44OTYsMTAsNCwxMHogTTI4LDE0SDRjLTEuMTA0LDAtMiwwLjg5Ni0yLDIgIHMwLjg5NiwyLDIsMmgyNGMxLjEwNCwwLDItMC44OTYsMi0yUzI5LjEwNCwxNCwyOCwxNHogTTI4LDIySDRjLTEuMTA0LDAtMiwwLjg5Ni0yLDJzMC44OTYsMiwyLDJoMjRjMS4xMDQsMCwyLTAuODk2LDItMiAgUzI5LjEwNCwyMiwyOCwyMnoiLz48L3N2Zz4=';
  my $closedMenuBase64 = 'data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiA/PjxzdmcgZmlsbD0ibm9uZSIgaGVpZ2h0PSIyNCIgc3Ryb2tlPSIjMDAwIiBzdHJva2UtbGluZWNhcD0icm91bmQiIHN0cm9rZS1saW5lam9pbj0icm91bmQiIHN0cm9rZS13aWR0aD0iMiIgdmlld0JveD0iMCAwIDI0IDI0IiB3aWR0aD0iMjQiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+PGxpbmUgeDE9IjE4IiB4Mj0iNiIgeTE9IjYiIHkyPSIxOCIvPjxsaW5lIHgxPSI2IiB4Mj0iMTgiIHkxPSI2IiB5Mj0iMTgiLz48L3N2Zz4=';

  $html .= "<html>\n<head><title>SneakerNet report</title>\n";

  # meta tags are in head
  $html.='<meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">'."\n";
  # https://purecss.io/start/#add-the-viewport-meta-element
  $html.='<meta name="viewport" content="width=device-width, initial-scale=1">'."\n";

  $html .= "<!-- Favicon image credit: shoe print by John Winowiecki from the Noun Project -->\n";
  $html .= "<link rel='icon' href='data:image/gif;base64,R0lGODlhDQAQAPUsAAAAAAEBAQICAgMDAw4ODhcXFx4eHi8vLzIyMjQ0NDo6OkZGRklJSU9PT1FRUVhYWFpaWmRkZGVlZWlpaXNzc4KCgoWFhYaGhpOTk5aWlp6enqCgoL6+vsLCwsfHx83NzdbW1t/f3+np6erq6uvr6+zs7O3t7fHx8fPz8/f39/v7+/7+/v///wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAAAAAAALAAAAAANABAAAAZMQJaQZaI8PsOkMDMAMFBKocoBABA6URZKUQVksqlFt5JlWbqWMshQ3ZSnVk+ZpXlMRuWRIxGBZkEFAAd4WSMIABEqcxwWIXOPkJFCQQA7'>\n";

  $html .= "<style>\n";
  $html .= "h1    {font-weight:bold; font-size:24px; margin:3px 0px;}\n";
  $html .= "h2    {font-weight:bold; font-size:16px; margin:3px 0px;}\n";
  $html .= "body  {font-size:12px; font-family:-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif, Apple Color Emoji, Segoe UI Emoji, Segoe UI Symbol;}\n";
  $html .= "table {border:1px solid black;}\n";
  $html .= "td    {border:1px solid #999999; margin:0px; word-wrap: break-all;}\n";
  $html .= "th    {border:1px solid #009900; margin:0px; word-wrap: break-all;}\n";
  $html .= "thead {font-weight:bold; background-color:#BBE;}\n";
  $html .= "tbody {color:black;} \n";
  $html .= "tfoot {font-size:12px; color:#DD3333;}\n";
  $html .= ".genericInfo {background-color:#EEEEEE; border: 1px solid #666666; margin:2px 0px; padding:1px;}\n";
  $html .= ".genericInfo p {margin-bottom:1em;}\n";
  $html .= ".version {font-size:12px; padding:1px; margin-top:10px; font-family:monospace; font-size:10px;}\n";
  $html .= ".pluginSplash{background-color:#FAFFFF;border:1px solid black; margin:6px 0px; padding:1px;}\n";
  # https://codeburst.io/how-to-make-a-collapsible-menu-using-only-css-a1cd805b1390
  $html .= ".pluginContent{display:none;}\n";
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
  $html .= "input:checked ~ .pluginContent{display:block;}\n";
  $html .= "input:checked ~ label {
    background-image: url(\"$closedMenuBase64\");
    transition: all ease 0.1s;
  }\n";
  $html .= "a { 
    color:black;
    text-decoration:none;
  }\n";
  $html .= "a:link { 
    color:red;
  }\n";
  $html .= "a:visited { 
    color:green;
  }\n";
  $html .= "a:hover { 
    color:hotpink;
  }\n";
  $html .= "a:hover { 
    color:blue;
  }\n";
  $html .= "</style>\n";

  $html .= "</head>\n";
  
  return $html;
}

sub htmlFooters{
  my($dir,$settings)=@_;

  my $realpath = realpath($dir);
  my $snPath   = $realpath . '/SneakerNet';

  my $html = "";

  $html .= "<p>For more information, see <a href='https://github.com/lskatz/SneakerNet'>https://github.com/lskatz/SneakerNet</a></p>\n";
  $html .= "<p>SneakerNet files can be found in <span style='font-weight:bold'>$snPath</span></p>\n";
  $html .= "</html>\n";

  return $html;
}

sub usage{
  print "Make an HTML report from SneakerNet results
  Usage: $0 MiSeq_run_dir
  --version
";
  exit(0);
}

