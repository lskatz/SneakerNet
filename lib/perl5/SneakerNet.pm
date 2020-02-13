package SneakerNet;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
#use File::Spec ();
use Config::Simple;
use Data::Dumper;
use Carp qw/croak confess carp/;

use FindBin qw/$Bin $Script $RealBin $RealScript/;

our @EXPORT_OK = qw(
  readConfig samplesheetInfo samplesheetInfo_tsv passfail
  command logmsg fullPathToExec version recordProperties readProperties
  exitOnSomeSneakernetOptions
);

our $VERSION = '0.8.1';

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
    croak $errStr;
  }
	return $fullpath;
}

# Takes a hash of executables => `way to check version`
# Some keys in the hash however are special and have an
# underscore at the beginning: _CITATION and _VERSION.
# Exits with 0 if any uses are invoked.
sub exitOnSomeSneakernetOptions{
  my($properties, $settings) = @_;

  if($$settings{version}){
    print $$properties{_VERSION}."\n";
    exit 0;
  }
  if($$settings{citation}){
    print $$properties{_CITATION}."\n";
    exit 0;
  }
  if($$settings{'check-dependencies'}){
    
    logmsg "$0: ".$$properties{_VERSION};
    my @exe = sort(
      grep {!/^_/}
      keys(%$properties)
    );

    # Print off all executable names before possible errors
    # down below.
    # Prints on stdout.
    for my $exe(@exe){
      print "$exe\n";
    }

    # Run through all execs but die if not present.
    # Prints on stderr
    for my $exe(@exe){
      my $path = fullPathToExec($exe);

      my $ver = 'UNKNOWN VERSION';
      if(my $vcmd = $$properties{$exe}){
        ($ver) = qx($vcmd);
        $ver or die "ERROR: could not determine version of '$exe' via '$vcmd'";
        chomp($ver);
        logmsg "$exe: $ver";
      }
    }

    exit(0);
  }

  return 0;
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

# Read the sneakernet-style spreadsheet
sub samplesheetInfo_tsv{
  my($samplesheet,$settings)=@_;

  # Get possible taxon rules.
  my $config = readConfig();
  
  my %sample;
  open(my $fh, "<", $samplesheet) or croak "ERROR: reading $samplesheet";
  while(<$fh>){
    chomp;
    my @F = split /\t/;
    my($sampleName,$rules,$fastq)=@F;
    $fastq ||= "";
    my @fastq = split(/;/, $fastq);
    $sample{$sampleName}={
      fastq => \@fastq,
      sample_id => $sampleName,
    };
    for my $rule(split(/;/, $rules)){
      my($key,$value)=split(/=/,$rule);
      $key=lc($key); # ensure lowercase keys
      my @values = split(/,/,$value);
      if(@values > 1){
        $sample{$sampleName}{$key} = \@values;
      } else {
        $sample{$sampleName}{$key} = $value;
      }
    }

    # Set up the taxon rules if possible
    $sample{$sampleName}{taxonRules}={};
    if(defined(my $taxon = $sample{$sampleName}{taxon})){
      my $possibleRules = $$config{obj}{"taxonProperties.conf"}->param(-block=>$taxon);
      if(defined($possibleRules)){
        $sample{$sampleName}{taxonRules}=$possibleRules;
      }
    }
  }
  close $fh;

  return \%sample;
}

sub samplesheetInfo{
  my($samplesheet,$settings)=@_;

  # If this is a tsv file, it is the simplified kind
  if($samplesheet=~/\.tsv$/){
    return samplesheetInfo_tsv(@_);
  }

  my $config=readConfig();

  my $section="";
  my @header=();
  my %sample;
  open(SAMPLE,$samplesheet) or croak "ERROR: could not open sample spreadsheet $samplesheet: $!";
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
      # Force the description to be empty string if:
      #   not defined
      #   Excel has '#N/A'
      if(!$F{description} || $F{description}=~/#N\/?A/i){
        $F{description}="";
      }

      for my $keyvalue(split(/;/,lc($F{description}))){
        my($key,$value)=split(/=/,$keyvalue);
        $value||="";
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
      croak "ERROR: could not find sample id for this line in the sample sheet: ".Dumper \%F if(!$F{sample_id});

      # What rules under taxonProperties.conf does this
      # genome mostly align with?
      my $alignedWith="";
      my %taxonProperties=%{ $$settings{obj}{"taxonProperties.conf"}->vars };
      my @taxa=$$settings{obj}{"taxonProperties.conf"}->get_block;
      #croak Dumper $$settings{obj}{"taxonProperties.conf"}->param(-block=>'Salmonella');
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
        if($F{taxonRules}{dest_subfolder}){
          logmsg "WARNING: dest_subfolder was set for UNKNOWN in taxonProperties.conf. This parameter is deprecated for UNKNOWN and will be ignored in the future. Please remove dest_subfolder under UNKNOWN in taxonProperties.conf and instead set catchall_subfolder in settings.conf.";
        }
        $F{taxonRules}{dest_subfolder} ||= $$settings{catchall_subfolder};
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
    if(!@possibleFastq){
      logmsg "WARNING: there is a sample $samplename but no files $samplename*.fastq.gz";
    }
    #for (@possibleFastq){
    #  $_ = File::Spec->abs2rel($_, dirname($samplesheet));
    #}
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
    croak $msg;
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
  open(my $passfailFh, $passfail) or croak "ERROR: could not read $passfail: $!\n  Please make sure that sn_passfail.pl is run before this script, but after the read metrics script.";
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
  return $VERSION;
}

# Record the plugin version and any other misc things
# into a run directory.
# Returns length of string that was written.
sub recordProperties{
  my($runDir,$writeHash, $settings)=@_;

  my $propertiesFile = "$runDir/SneakerNet/properties.txt";
  my $writeString="";
  if(!-e $propertiesFile || (stat($propertiesFile))[7] == 0){
    $writeString.=join("\t", qw(plugin key value))."\n";
  }
  for my $key(keys(%$writeHash)){
    if(!defined($$writeHash{$key})){
      carp "WARNING: in SneakerNet::recordProperties(), key '$key' was not defined";
      next;
    }
    $writeString.=join("\t",basename($0), $key, $$writeHash{$key})."\n";
  }

  open(my $fh, ">>", $propertiesFile) or croak "ERROR writing to $propertiesFile: $!";
  print $fh $writeString;
  close $fh;
  
  return length($writeString);
}

# Read properties, the opposite of recordProperties().
# Returns a properties hash of hash, where the primary
# key is the plugin, and each plugin has a hash.
# Each plugin should have a "version" key/value.
# E.g., $property{"guessTaxon.pl"}{version} = 1.0
sub readProperties{
  my($runDir, $settings) = @_;
  my %prop = ();
  my $propertiesFile = "$runDir/SneakerNet/properties.txt";
  open(my $fh, '<', $propertiesFile) or croak "ERROR reading $propertiesFile: $!";
  my $header = <$fh>;
  while(my $line = <$fh>){
    chomp($line);
    my($plugin, $key, $value) = split(/\t/, $line);
    $prop{$plugin}{$key} = $value;
  }

  return \%prop;
}


1;

