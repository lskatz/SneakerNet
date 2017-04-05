package SneakerNet;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use Config::Simple;
use Data::Dumper;
use Carp qw/croak confess/;

use FindBin qw/$Bin $Script $RealBin $RealScript/;

our @EXPORT_OK = qw(
  readConfig samplesheetInfo passfail
  command logmsg fullPathToExec version
);


my $thisdir=dirname($INC{'SneakerNet.pm'});

sub logmsg{print STDERR "$0: @_\n";}

# If argument is an executable in the current path, returns 
# the full path to it, otherwise dies.
sub fullPathToExec($;$) {
	my ($executable,$settings) = @_;
	my $fullpath="";
	for ('',split(/:/, $ENV{PATH})) {
    my $path=$_."/".$executable;
		if (-x $path && -f $path) { $fullpath = File::Spec->rel2abs($path); last; }
    if (-x $path && -l $path && -f readlink($path)){ $fullpath = File::Spec->rel2abs(readlink($path)); last; }
	}
  if(! -x $fullpath){
	  my $errStr="Error finding full path to executable ($executable)";
    die $errStr;
  }
	return $fullpath;
}

sub readConfig{
  my $settings={};

  my @file=glob("$thisdir/../../config/*.conf");
  for my $file(@file){
    my $cfg = new Config::Simple();
    if(!$cfg->read($file)){
      logmsg "WARNING: could not read $file: ".$cfg->error;
      next;
    }
    my %vars= $cfg->vars();
    $$settings{$_}=$vars{$_} for(keys(%vars));
    $$settings{obj}{basename($file)}=$cfg; # save the obj too
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
      $F{route}||=[]; # force route to be an array
      $F{description}||="";
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

      # The HiSeq seems to use 'sampleid' instead of 'sample_id'
      if(!$F{sample_id}){
        $F{sample_id}=$F{sampleid};
      }
      die "ERROR: could not find sample id for this line in the sample sheet: ".Dumper \%F if(!$F{sample_id});

      # What rules under taxonProperties.conf does this
      # genome mostly align with?
      my $alignedWith="";
      my %taxonProperties=%{ $$settings{obj}{"taxonProperties.conf"}->vars };
      my @taxa=$$settings{obj}{"taxonProperties.conf"}->get_block;
      #die Dumper $$settings{obj}{"taxonProperties.conf"}->param(-block=>'Salmonella');
      for my $taxon(@taxa){
        my $taxonRegex=$$settings{obj}{"taxonProperties.conf"}->param("$taxon.regex");
        
        if(
              $F{sample_id} =~ /$taxonRegex/
           || ($F{species} && $F{species}=~/$taxon/i)
          ){

          $F{taxonRules}=$$settings{obj}{"taxonProperties.conf"}->param(-block=>$taxon);
          $F{taxonRules}{taxon}=$taxon;
          $F{species}=$taxon;
          last;
        }
      }
      # set it to the unknown if it's still not known
      if(!$F{species}){
        $F{species}="UNKNOWN";
        $F{taxonRules}=$$settings{obj}{"taxonProperties.conf"}->param(-block=>$F{species});
        $F{taxonRules}{taxon}=$F{species};
      }

      $sample{$F{sample_id}}=\%F;
      
    }
  }
  close SAMPLE;

  if(keys(%sample) == 0){
    logmsg "WARNING: there were zero samples found in the sample sheet. Is there a section labeled [data]?\n  in $samplesheet";
  }

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
  my $stdout=`$command`;
  if($?){
    my $msg="ERROR running command\n  $command";
    confess $msg;
  }

  return $stdout;
}

# Which files failed?
sub passfail{
  my($dir,$settings)=@_;

  # Which files should be skipped according to Q/C?
  # Read the passfail file which should have Sample as a
  # header and then the rest of the headers are pass
  # or fail values.
  my $passfail="$dir/SneakerNet/forEmail/passfail.tsv";
  my %failure;
  open(my $passfailFh, $passfail) or die "ERROR: could not read $passfail: $!\n  Please make sure that sn_passfail.pl is run before this script, but after the read metrics script.";
  my $header=<$passfailFh>;
  chomp($header);
  my @header=split(/\t/,$header);
  while(<$passfailFh>){
    next if(/^#/);
    chomp;
    my @F=split(/\t/,$_);
    my %F;
    @F{@header}=@F;

    # Remove the sample header so that all values of 
    # %failure have to do with pass/fail
    my $sample=$F{Sample};
    delete($F{Sample});
    $failure{$sample}=\%F;
  }
  close $passfailFh;
  
  return \%failure;
}

# Return the version of SneakerNet
sub version{

  my $codeRepoVer="-1";
  my $configVer="-1";

  my $cfg = new Config::Simple();
  if(!$cfg->read("$thisdir/../../config.bak/settings.conf")){
    logmsg "WARNING: could not read $thisdir/../../config.bak/settings.conf: ".$cfg->error;
  }
  $codeRepoVer=$cfg->param("version");

  # See if the code's version matches the custom version
  my %settings=%{ readConfig() };
  $configVer=$settings{version} if($settings{version});

  if($configVer ne $codeRepoVer){
    logmsg "WARNING: the codebase version is reported differently than the configuration. Please review the config folder to update any new options and to update the version number.";
    logmsg "The current code repository version is $codeRepoVer.  The custom version is $configVer";
  }

  return $codeRepoVer;
}

1;

