package SneakerNet;
use strict;
use warnings;
use Exporter qw(import);
use File::Basename qw/fileparse basename dirname/;
use File::Copy qw/mv/;
use File::Temp qw/tempdir/;
#use File::Spec ();
use Config::Simple;
use Data::Dumper;
use Carp qw/croak confess carp/;
use Cwd qw/realpath/;
use List::Util qw/sum min max/;
use version 0.77;

use FindBin qw/$Bin $Script $RealBin $RealScript/;

our @EXPORT_OK = qw(
  readConfig samplesheetInfo samplesheetInfo_tsv passfail
  readTsv readKrakenDir
  command logmsg fullPathToExec version recordProperties readProperties
  exitOnSomeSneakernetOptions
  readMetrics
  @rankOrder %rankOrder $VERSION
);

=pod

=head1 NAME

SneakerNet is a QA/QC pipeline for a MiSeq/HiSeq/Ion Torrent/assembly-only run

=head1 SYNOPSIS

TODO

=cut

our $VERSION  = version->declare('0.26.0');
our %rankName = (S=>'species', G=>'genus', F=>'family', O=>'order', C=>'class', P=>'phylum', K=>'kingdom', D=>'domain', U=>'unclassified');
our @rankOrder= qw(S G F O C P K D U);
our %rankOrder= (S=>0, G=>1, F=>2, O=>3, C=>4, P=>5, K=>6, D=>7, U=>8);

my $thisdir=dirname($INC{'SneakerNet.pm'});

# Tools and how to detect their versions
# The ommands will generate the semver on a single line
# somewhere.
# We'd have to run the semver regex to really know the semver however.
our %toolVersion = (
  # GNU/Linux core
  cat         => 'cat --version | head -n 1',
  cp          => 'cp --version | head -n 1',
  "mv"        => 'mv --version | head -n 1',
  rm          => 'rm --version | head -n 1',
  sendmail    => 'apt-cache show sendmail 2>/dev/null | grep Version || rpm -qi sendmail 2>/dev/null | grep Version',
  uuencode    => 'uuencode --version | grep uuencode',
  sort        => 'sort --version | head -n 1',
  head        => 'head --version | head -n 1',
  uniq        => 'uniq --version | head -n 1',
  touch       => 'touch --version | head -n 1',
  zip         => 'zip --version | grep "This is Zip"',
  'python3'   => 'python3 --version',
  perl        => 'perl -v',
  wget        => 'wget --version | head -n 1',
  rsync       => 'rsync --version | grep version',
  ssh         => 'ssh -V 2>&1',

  # Misc binfie tools
  fastqc        => "fastqc --version",
  shovill       => 'shovill --version',
  prodigal      => "prodigal -v 2>&1 | grep -i '^Prodigal V'",
  'blastn'      => 'blastn -version | head -n 1',
  mlst          => 'mlst --version',
  'chewBBACA.py'=> 'chewBBACA.py -h 2>&1 | grep -m 1 -i version',
  samtools      => 'samtools 2>&1 | grep Version:',
  bcftools      => 'bcftools 2>&1 | grep Version:',
  'trimmomatic' => 'trimmomatic -version 2>&1 | grep -v _JAVA',
  seqtk         => 'seqtk 2>&1 | grep -m 1 Version:',
  # note: bgzip older versions are buggy for detecting
  # version and so this is a hack
  bgzip         => 'which bgzip && echo "0.0.1"',
  # buggy bgzip => 'bgzip --version | head -n1',
  tabix         => 'tabix --version | head -n1',
  bowtie2       => 'bowtie2 --version | grep -m 1 version',
  'bowtie2-build' => 'bowtie2-build --version | grep -m 1 version',
  staramr       => 'staramr --version 2>/dev/null',
  colorid       => 'colorid --version',
  kma           => 'kma -v',
  multiqc       => 'multiqc --version',
  quast         => 'quast --version',
  'SalmID.py'   => 'SalmID.py --version',

  # CGP doesn't have a simple version check
  # and so just have something idiotic
  'run_assembly_metrics.pl'     => "which run_assembly_metrics.pl >/dev/null 2>&1 && echo 0.0.1",
  'run_assembly_readMetrics.pl' => "which run_assembly_readMetrics.pl >/dev/null 2>&1 && echo 0.0.1",
  'run_prediction_metrics.pl'  => "which run_prediction_metrics.pl >/dev/null 2>&1 && echo 0.0.1",

  # shovill requires a ton of things:
  'seqtk'       => 'seqtk 2>&1 | grep Version',
  'pigz'        => 'pigz --version 2>&1',
  'mash'        => 'mash --version 2>&1',
  'pigz'        => 'pigz --version 2>&1',
  'trimmomatic' => 'trimmomatic -version 2>&1 | grep -v _JAVA',
  'lighter'     => 'lighter -v 2>&1',
  'flash'       => 'flash --version 2>&1 | grep FLASH',
  'spades.py'   => 'spades.py  --version 2>&1',
  'skesa'       => 'skesa --version 2>&1 | grep SKESA',
  'gfa_connector' => 'gfa_connector --version 2>/dev/null | grep gfa_connector',
  'bwa'         => 'bwa 2>&1 | grep Version:',
  'samtools'    => 'samtools 2>&1 | grep Version:',
  'samclip'     => 'samclip --version 2>&1',
  'java'        => 'java -version 2>&1 | grep version',
  'pilon'       => 'pilon --version 2>&1 | grep -v _JAVA',

  # not using megahit or velvet in this instance of shovill
  # but why not have these ready just in case
  'megahit'     => 'megahit --version 2>&1',
  'megahit_toolkit' => 'megahit_toolkit dumpversion 2>&1',
  'velveth'     => 'velveth 2>&1 | grep Version',
  'velvetg'     => 'velvetg 2>&1 | grep Version',

  # krona
  'ktImportKrona'     => 'ktImportKrona | grep "/" | grep -P -m 1 -o "KronaTools .*ktImportKrona"',
  'ktImportText'      => 'ktImportText | grep "/" | grep -P -m 1 -o "KronaTools .*ktImportText"',

  # kraken
  kraken              => 'kraken --version | grep -m 1 version',
  'kraken-translate'  => 'kraken-translate --version | grep -m 1 version',
  'kraken-report'     => 'kraken-report --version | grep -m 1 version',
);

=pod

=head2 PROPERTIES

=head3 $VERSION

The version of SneakerNet

=head3 %rankName

Taxonomy abbreviations, e.g., S=>species

=head3 @rankOrder

Array of rank abbreviations, sorted: S G F O C P K D U

=head3 %rankOrder

Hash of rank abbreviations and their order, e.g., S=>0, G=>1

Useful for comparing ranks
  
  use SneakerNet qw/$rankOrder/;
  $rankOrder{S} <=> $rankOrder{G};

=head2 METHODS

=head3 logmsg

This function prints a string to STDERR in the format of

    script: $arg1 $arg2 ...\n;

=cut

sub logmsg{print STDERR "$0: @_\n";}

=pod

=head3 fullPathToExec

fullPathToExec is a function that returns the full path
to an executable.
Croaks if the path does not exist.

Arguments:  executable (string)
returns:    absolute path (string)

=cut

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

=pod

=head3 versionOf

Find the specific version of an executable in the environment.
However, to find the version of a SneakerNet plugin, use the
plugin executable with `--version`.
Uses the `%toolVersion` hash that is global for this library.

Arguments: 

  $exe      => the executable, e.g., `blastn`.
  $settings => a hash of yet undetermined settings

Returns: semver-formatted version string or empty string if not found

=cut

sub versionOf{
  my($exe, $settings) = @_;
  
  # Suggested regex on https://semver.org/
  # Except I removed the anchors ^ and $
  # and I added a capture group for the whole thing which
  # should help avoid accidentally capturing only one
  # piece of the semver or another.
  my $semverRegex = qr/((?P<major>0|[1-9]\d*)\.(?P<minor>0|[1-9]\d*)\.(?P<patch>0|[1-9]\d*)(?:-(?P<prerelease>(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+(?P<buildmetadata>[0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?)/;

  my $cmd = $toolVersion{$exe}
    or die "ERROR: I do not know how to find the version of $exe. To add a new tool, go to \%toolVersion in SneakerNet.pm and add it.";
  my $version = `$cmd`;
  chomp($version);
  
  my $semver = "";

  if($version =~ $semverRegex){
    $semver = $1;
  } else{
    # If we see a major/minor number, 
    # try to add a ".0" for the patch number just in case that helps
    $version =~ s/(\d+\.\d+)/$1.0/;
    if($version =~ $semverRegex){
      $semver = $1;
    }
  }

  return $semver;
}

=pod

=head3 versionsOf

Find the specific versions of groups of executables.
Uses similarly named subroutine `versionOf` internally.

Arguments: 

  $group    => the executables group which could be one of the following
  $settings => a hash of yet undetermined settings

Returns: list of semver-formatted version strings

Named groups:

=over

=item all - get all available software versions

=back

=cut

sub versionsOf{
  my($group, $settings) = @_;
  my @semver;

  # Set to lowercase to make str comparisons easy
  my $lcGroup = lc($group);

  if($lcGroup eq 'all'){
    for my $exe(sort(keys(%toolVersion))){
      my $tmp = versionOf($exe, $settings);
      push(@semver, "$exe $tmp");
    }
  }

  # See if this is a workflow
  else {
    my $snSettings = readConfig();
    my $workflows = $$snSettings{obj}{'plugins.conf'}->param(-block=>'plugins');

    my $plugins = $$workflows{$lcGroup};
    if(defined($plugins)){
      for my $p(@$plugins){
        my @version = `$p --check-dependencies 2>&1 >/dev/null`;
        for my $v(@version){
          chomp($v);
          if($v =~ /([\w\.]+):\s+(\w+)\s+(.+)/){
            my $plugin = $1;
            my $key = join("-", $2, "version");
            my $value = $3;
            push(@semver, "$plugin $key $value");
          }
        }
      }
    } else {
      carp "ERROR: could not determine what to get versions of with keyword $lcGroup (did you try to use an undefined workflow name?";
    }
  }


  return \@semver;
}

=pod

=head3 exitOnSomeSneakernetOptions

exitOnSomeSneakernetOptions Takes a hash in the form of {executables => `way to check version`}
Some keys in the hash however are special and have an
underscore at the beginning.

Arguments: 

  $properties (hash ref)
    executable         => string (command that lets you check a version. Can be in the format of 'executable ...' where anything after a space is ignored for determining the viability of the executable)
    _CITATION          => bool   (prints $$settings{citation} and exits 0)
    _VERSION           => bool   (prints $$settings{version} and exits 0)
  $settings   (hash ref)
    version            => bool
    citation           => string
    check-dependencies => bool

Returns: 0

=cut

sub exitOnSomeSneakernetOptions{
  my($properties, $settings) = @_;
  $$properties{exe} //= [];

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
    my @exe = sort(@{ $$properties{exe} });

    # Print off all executable names before possible errors
    # down below.
    # Prints on stdout.
    for my $exe(@exe){
      print "$exe\n";
    }

    # Run through all execs but die if not present.
    # Prints on stderr
    my $numNotFound = 0;
    for(my $i=0;$i<@exe;$i++){
      my $desc = $exe[$i]; # avoid changing by ref
      my $exe  = $desc;
         $exe =~ s/\s+.*//; # remove anything after a whitespace
      my $path = eval{fullPathToExec($exe);};
      if(!defined($path) || !-e $path){
        logmsg "ERROR: could not find path to $exe";
        $numNotFound++;
        next;
      }

      my $ver = versionOf($exe, $settings);
      if($ver){
        logmsg "$exe $ver";
      } else {
        logmsg "ERROR: I could not find the version of $exe";
        $numNotFound++;
      }
      #my $ver = 'UNKNOWN VERSION';
      #if(my $vcmd = $toolVersion{$exe}){
      #  ($ver) = qx($vcmd);
      #  if(!$ver){
      #    logmsg "ERROR: could not determine version of '$desc' via '$vcmd'";
      #    $numNotFound++;
      #  } else {
      #    chomp($ver);
      #    logmsg "$exe: $ver";
      #  }
      #}

    }
    if($numNotFound > 0){
      croak "$numNotFound dependencies were not found.";
    }

    exit(0);
  }

  return 0;
}

=pod

=head3 readConfig

readConfig reads all files ending with .conf in the installation folder
under config/ (ie, SneakerNet/config/*.conf) using L<Config::Simple>.

Arguments: none

Returns: 
  
  hash reference with
    key=>value from all key/values in conf files
    Config::Simple objects under the obj key. Each key under obj is a basename of the config file.
      For example, this key/value will exist: $$hash{obj}{settings}{KRAKEN_DEFAULT_DB} = "/path/to/kraken/db"
      Example of comma-separated list in a conf file: $$hash{obj}{settings}{happiness_range} = ["&#128515;", "&#128556;", "&#128561;"]

=cut

sub readConfig{
  my $settings={};

  #my %seenBasename = ();

  my @file=glob("$thisdir/../../config/*.conf $thisdir/../../config.bak/*.conf");
  for my $file(@file){
    # Look at the filename. If we've seen this filename,
    # then skip. This will happen if the file exists
    # in config/ and in config.bak/
    #next if($seenBasename{basename($file)}++);

    my $cfg = new Config::Simple();
    if(!$cfg->read($file)){
      logmsg "WARNING: could not read $file: ".$cfg->error;
      next;
    }
    my %vars= $cfg->vars();
    for my $key(keys(%vars)){
      if(dirname($file) =~ /config.bak/){
        if(!defined($$settings{$key})){
          logmsg "WARNING: key $key was expected but was not found in your local config. I will use the value found in $file. Please edit in $thisdir/../../config/";
        }
      }

      $$settings{$key} //= $vars{$key};
    }
    $$settings{obj}{basename($file)} //= $cfg; # save the obj too
  }
  return $settings;
}

=pod

=head3 readKrakenDir

Read a Kraken directory from sn_kraken.pl

Arguments:

  directory:  Directory of kraken results, e.g., "dir/SneakerNet/kraken/sample1"
  settings:   Reference of a hash with the following keys:
    minpercent:  minimum percentage to consider for inclusion (float 0 to 100)

Returns:

  Results hash with keys corresponding to @SneakerNet::rankOrder
  Each key is an array of hashes, e.g.,
   C => [
          {
            'percent' => '75.31',
            'taxname' => 'Gammaproteobacteria',
            'rank' => 'C',
            'classifiedReads' => '68077',
            'taxid' => '1236',
            'specificClassifiedReads' => '40'
          },
          {
            'percent' => '24.49',
            'taxname' => 'Betaproteobacteria',
            'rank' => 'C',
            'classifiedReads' => '22141',
            'taxid' => '28216',
            'specificClassifiedReads' => '27'
          }
        ];

  On failure, returns empty hash.

=cut

sub readKrakenDir{
  my($dir, $settings) = @_;

  my $minPercent = $$settings{minpercent};
  if(!defined($minPercent)){
    $minPercent = 1e-6;
    logmsg "WARNING: minpercent was not set for readKrakenDir(). Setting to $minPercent";
  }
  
  my %bestGuess;
  my $taxfile = "$dir/kraken.report";
  if(! -e $taxfile){
    logmsg "WARNING: kraken.report does not exist at $taxfile - returning empty hash.";
    return {};
  }
  my @header = qw(percent classifiedReads specificClassifiedReads rank taxid taxname);
  open(my $fh, '<', $taxfile) or die "ERROR: could not open $taxfile for reading: $!";
  while(my $taxline = <$fh>){
    # trim whitespace
    $taxline =~ s/^\s+|\s+$//g;

    # Split fields on tab, but also remove whitespace on fields
    my @F = split(/\s*\t\s*/, $taxline);

    # Name the fields
    my %field;
    @field{@header} = @F;

    # We only care about named ranks like phylum, family, genus, or species
    next if($field{rank} eq '-');

    # Filter for low percentage counts
    next if($field{percent} < $minPercent);

    push(@{ $bestGuess{$field{rank}} }, \%field);

  }
  close $fh;

  # Taxonomy ranks in order
  my @sortedRank = qw(S G F O C P K D U);
  # Index the ranks so that we can sort, e.g., $rank{S} <=> $rank{F}
  my %rank;
  for(my $i=0;$i<@sortedRank;$i++){
    $rank{$sortedRank[$i]} = $i;
  }

  # sort %bestGuess within each rank, so that the best guess
  # is always element zero for each rank.
  my %raw;
  while(my($rank, $taxa) = each(%bestGuess)){
    my @sortedTaxon = sort{$$b{percent} <=> $$a{percent} } @$taxa;
    $raw{$rank} = \@sortedTaxon;
  }

  # Ensure that each rank is filled in, even if it has to
  # be a zero percentage.
  for my $rank(@sortedRank){
    if(ref($raw{$rank}) ne 'ARRAY'){
      $raw{$rank} = [ {percent=>0, taxname=>"UNKNOWN", rank=>$rank, classifiedReads=>0, taxid=>-1, specificClassifiedReads=>0} ];
    }
  }

  return \%raw;
}

=head3 readTsv

Reads a generic TSV (tab separated values) file

Arguments:

  filename:  path to tsv-formatted file
  settings:  hash of settings with keys:
    headers:   Whether or not the first line is headers (bool, default: true)
    keyIndex:  Which column is the key for each row (int, default: 0)

  Returns:   reference to a hash of values

=cut

sub readTsv{
  my($filename, $settings) = @_;

  my ($package, $script, $line, $subroutine, $hasargs,$wantarray, $evaltext, $is_require, $hints, $bitmask, $hinthash) = caller(1);

  my %TSV;

  # If the filename can't be read or if it has zero bytes, return empty hash
  if(!-r $filename){
    logmsg "$subroutine (line $line): WARNING: could not read from $filename";
    return \%TSV;
  }
  if(!-s $filename){
    logmsg "$subroutine (line $line): WARNING: trying to read from zero byte file $filename";
    return \%TSV;
  }

  my $keyIndex = $$settings{keyIndex} || 0;
  my $headers  = 1;
  if(defined($$settings{headers}) && $$settings{headers}==0){
    $headers = 0;
  }
  
  open(my $fh, '<', $filename) or croak("ERROR: could not read from TSV file $filename: $!");
  my @header;
  if($headers){
    my $header = <$fh>;
    chomp($header);
    @header = split(/\t/, $header);
  }
  while(<$fh>){
    chomp;
    my @F = split /\t/;

    # If the header isn't set, then it is a range of ints
    if(!$headers){
      @header = (0 .. scalar(@F)-1);
    }

    # Map keys=>values
    my %F;
    @F{@header} = @F;

    my $rowKey = $F[$keyIndex];
    $TSV{$rowKey} = \%F;
  }
  close $fh;

  return \%TSV;
}


=head3 samplesheetInfo_tsv

Reads a samples.tsv file in a SneakerNet run directory.

Arguments:

  Samples file: path/to/samples.tsv (string)

Returns:

  $samples: information on each sample (hash ref)
    sampleName  => {
      taxonProperties=>{
        coverage => 20,
        ...
      },
      taxon => "Listeria",
      fastq => [
        path/to/R1,
        path/to/R2,
      ],
      ...
    }, 
    sampleName2 => { ... },

=cut

sub samplesheetInfo_tsv{
  my($samplesheet,$settings)=@_;

  my $runDir = realpath(dirname($samplesheet));

  # Get possible taxon rules.
  my $config = readConfig();

  # If species have been detected already, load them up.
  # The kraken folder will be the way we do that for now.
  my %speciesSuggestion;
  my $krakenReport = "$runDir/SneakerNet/forEmail/kraken.tsv";
  if(-e $krakenReport){
    open(my $fh, '<', $krakenReport) or die "ERROR: could not read from $krakenReport: $!";
    my $header = <$fh>;
    chomp($header);
    my @header = split(/\t/, $header);
    while(<$fh>){
      next if(/^#/);
      chomp;
      my %F;
      my @F = split /\t/;
      @F{@header} = @F;
      my @suggestions = ($F{BEST_GUESS}, (split(/\s+/, $F{BEST_GUESS}))[0]);
      $speciesSuggestion{$F{NAME}} = \@suggestions;
    }
    close $fh;
  }
  # TODO species suggestion from other sources???
  #print Dumper \%speciesSuggestion;

  my %sample;
  open(my $fh, "<", $samplesheet) or croak "ERROR: reading $samplesheet";
  while(<$fh>){
    chomp;
    next if(/^\s*$/);
    my @F = split /\t/;
    my($sampleName,$rules,$fastq)=@F;
    $fastq ||= "";
    my @fastq = split(/;/, $fastq);
    # get the abs path to each fastq file
    for(@fastq){
      $_ = File::Spec->rel2abs($_, $runDir);
    }
    $sample{$sampleName}={
      fastq => \@fastq,
      sample_id => $sampleName,
    };

    # Find out whether there is an assembly for this sample
    my $asmDir = "$runDir/SneakerNet/assemblies/$sampleName";
    $sample{$sampleName}{asm} = "";
    if(-d $asmDir){
      my @asm = glob("$asmDir/*.fasta");
      if(@asm){
        $sample{$sampleName}{asm} = $asm[0];
      }
    }

    my @rule = grep{$_ ne '.'} split(/;/, $rules);
    for my $rule(@rule){
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
    # Get a few options for the taxon: explicit from the samplesheet or guessed
    my $explicitTaxon = $sample{$sampleName}{taxon};
    my $guessedTaxons = $speciesSuggestion{$sampleName} || [];
    my @taxonGuess = grep {defined($_)} ($explicitTaxon, @$guessedTaxons);
    @taxonGuess    = ((grep {!/^unknown$/i} @taxonGuess), 'UNKNOWN');
    # Loop through the different taxa guesses, starting with
    # what was explicitly given
    for my $taxon(@taxonGuess){
      if(defined($taxon) && $taxon ne ""){
        my $possibleRules = $$config{obj}{"taxonProperties.conf"}->param(-block=>$taxon);
        # This is the right taxon if we have rules associated with it
        if(defined($possibleRules) && scalar(keys(%$possibleRules))>1){
          $sample{$sampleName}{taxonRules}=$possibleRules;
          $sample{$sampleName}{taxon} = $taxon;
          # If we give it taxon rules, then we're done going through different taxon guesses
          last;
        }
      }
    }

    # Set the file path to the reference fasta if it exists
    my $ref_id = $sample{$sampleName}{taxonRules}{reference_fasta_id} || [];
    $ref_id = [$ref_id] if(ref($ref_id) ne 'ARRAY');
    my $ref_name = join("", @$ref_id);
    my $ref_commas=join(",", @$ref_id);
    if(defined($ref_id)){
      my $dir = realpath($RealBin."/../db/fasta");
      mkdir($dir);
      my $ref = "$dir/$ref_name.fasta";
      my $gbk = "$dir/$ref_name.gbk";

      if(!-e $ref){
        logmsg "Did not find @$ref_id in the SneakerNet installation. Downloading @$ref_id into $gbk and $ref";
        logmsg "  Wget log can be found at $dir/*.log";
        mkdir dirname($dir);
        mkdir $dir;
        my $wgetxopts = "";
        if($ENV{NCBI_API_KEY}){
          $wgetxopts .= "&ncbi_api_key=$ENV{NCBI_API_KEY}";
        }
        command("wget 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$ref_commas&rettype=gbwithparts&format=genbank$wgetxopts' -O $gbk 2> $gbk.log ");

        # The file is downloaded into a .tmp file and then moved,
        # so that it is clear when it is 100% downloaded.
        command("wget 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=$ref_commas&rettype=fasta$wgetxopts' -O $ref.tmp 2> $ref.log ");
        mv("$ref.tmp", $ref);
      }

      $sample{$sampleName}{taxonRules}{reference_gbk}   = $gbk;
      $sample{$sampleName}{taxonRules}{reference_fasta} = $ref;
    }

    # Get the primers bed file if it doesn't exist
    my $primers = $sample{$sampleName}{taxonRules}{primers_bed_url};
    if(defined($primers)){
      my $dir = realpath($RealBin."/../db/fasta");
      my $bed = "$dir/".basename($primers);
      if(!-e $bed){
        logmsg "Did not find local bed and so I am downloading $primers to $bed";
        mkdir dirname($dir);
        mkdir $dir;
        command("wget '$primers' -O $bed");
      }
      $sample{$sampleName}{taxonRules}{primers_bed} = $bed;
    }
    $sample{$sampleName}{taxonRules}{primers_bed} //= "";
    $sample{$sampleName}{taxonRules}{primers_bed} =~ s|__BASE__|$thisdir/../..|g;
  }
  close $fh;

  return \%sample;
}

=head3 samplesheetInfo

This should be rarely used outside of F<sn_parseSampleSheet.pl>. Parses an
Illumina-style SampleSheet.csv file. If the file ends in .tsv, runs
L</samplesheetInfo_tsv>.

Arguments: path/to/SampleSheet.csv

Returns: $sample (hash ref). For more information on $sample: L</samplesheetInfo_tsv>.

=cut

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

    # If we see this syntax, then we change $section
    if(/^\[(\w+)\]/){  # [sectionname]
      $section=lc($1);
      # store into @header the field names
      my $header=<SAMPLE>;
      $header=~s/^\s+|\s+$//g; # trim whitespace
      @header=split(/,/,lc($header));
      next;
    }

    # This is where the data are stored for each sample
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

      # Get the key/value pairs in that description field.
      # They are separated by semicolon, and key/values by equals.
      for my $keyvalue(split(/;/,lc($F{description}))){
        my($key,$value)=split(/=/,$keyvalue);
        $value||="";
        $key=~s/^\s+|\s+$//g;      #whitespace trim
        $value=~s/^\s+|\s+$//g;    #whitespace trim
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
      # This usually just skips blank lines
      if(!$F{sample_id}){
        carp "SKIP: could not find sample id for this line in the sample sheet: ".Dumper \%F;
        carp "If you think this SKIP was in error, please fill in the sample_id under SampleSheet.csv and run this script again with --force.";
        next;
      }

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

=head3 command

Runs a system command. STDERR is not captured and will leak through to be printed.
If a warning is produced with $?, will croak with the error message.

Arguments:

  $command:  a command that will be used in L</system> (string)
  $settings: controls for how the command will be run (hash ref)
    debug  => bool (repeat the command back to the terminal with L</logmsg>)

Returns:

  $stdout:   standard output from command (string)

=cut

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

=head3 passfail

Reads F<passfail.tsv> in the SneakerNet results directory.

Arguments:  SneakerNet run directory path

Returns:

  $failure (hash ref)
    # 0 for pass, 1 for fail, -1 for unknown
    sample1 => {
      coverage => -1,
      quality  => 1,
      kraken   => 0,
    },
    sample2 => {
      coverage => 0,
      quality  => 0,
      kraken   => 0,
    },
    sample3 => {...},

=cut
  
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

=head3 version

Returns the version of SneakerNet, encoded under $SneakerNet::Version

=over 

=item $SneakerNet::VERSION

$SneakerNet::VERSION is the current version of SneakerNet

=back

=cut

sub version{
  return $VERSION;
}

=head3 recordProperties

Record the plugin version and any other misc things
into a run directory.
Returns length of string that was written.
Writes to F<properties.txt> as a rudimentary local database

Arguments:

  $runDir:    The path to the SneakerNet run directory (string)
  $writeHash: Hash of properties to record. Any property name is viable. Common properties:
    version   =>  string (version number)
    table     =>  string (path to a results file in tsv format)
    date      =>  string (date of when results finished in YYYY-MM-DD format
    time      =>  string (time of when results finished in HH:MM:SS format
    exe       =>  list of executables that can be defined in `versionOf()`

Returns:  number of characters written to F<properties.txt>

=cut

sub recordProperties{
  my($runDir,$writeHash, $settings)=@_;

  my $propertiesFile = "$runDir/SneakerNet/properties.txt";
  if(!-d "$runDir/SneakerNet"){
    mkdir "$runDir/SneakerNet";
  }

  my $writeString="";
  if(!-e $propertiesFile || (stat($propertiesFile))[7] == 0){
    $writeString.=join("\t", qw(plugin key value))."\n";
  }
  for my $key(keys(%$writeHash)){
    if(!defined($$writeHash{$key})){
      carp "WARNING: in SneakerNet::recordProperties(), the value for key '$key' was not defined";
      next;
    }
    if($key eq 'exe'){
      if(ref($$writeHash{$key}) ne 'ARRAY'){
        confess "ERROR: key $key does not contain an array";
      }
      for my $exe(sort @{ $$writeHash{$key} }){
        $writeString .= join("\t", basename($0),
          "$exe-version", versionOf($exe, $settings))."\n";
      }
    }
    # If the key isn't special like 'exe' then it's a string and just
    # treat it normally
    else {
      $writeString.=join("\t",basename($0), $key, $$writeHash{$key})."\n";
    }
  }

  open(my $fh, ">>", $propertiesFile) or croak "ERROR writing to $propertiesFile: $!";
  print $fh $writeString;
  close $fh;
  
  return length($writeString);
}

=head3 readProperties

Read properties, the opposite of L</recordProperties>.
Returns a properties hash of hash, where the primary
key is the plugin, and each plugin has a hash.
Each plugin should have a "version" key/value.

Arguments:

  $runDir:    The path to the SneakerNet run directory (string)

Returns

  $properties:  hash ref.  E.g.,
    $property{"guessTaxon.pl"}{version} = "1.0"

=cut

sub readProperties{
  my($runDir, $settings) = @_;
  my %prop = ();
  my $propertiesFile = "$runDir/SneakerNet/properties.txt";
  open(my $fh, "tac $propertiesFile | ") or croak "ERROR reading $propertiesFile: $!";
  #my $header = <$fh>;
  while(my $line = <$fh>){
    next if($line =~ /^\s*$/); # skip blank lines
    chomp($line);
    my($plugin, $key, $value) = split(/\t/, $line);
    next if($plugin =~ /^plugin$/i); # skip the header

    next if(defined $prop{$plugin}{$key});
    if($key =~ /table/i){
      my $path = $value;
      if(!-f $path){
        $path = "$runDir/SneakerNet/forEmail/".basename($value);
      }
      if(!-f $path){
        die "ERROR ($plugin): value for table was given as $value but it was not found";
      }
      $path = File::Spec->rel2abs($path);

      $value = realpath($path);
      logmsg $value;
    }
    $prop{$plugin}{$key} = $value;
  }

  return \%prop;
}

=head3 readMetrics(\@fastq, $expectedGenomeSize, $settings)

Read metrics.
Replacing CG-Pipeline's run_assembly_readMetrics.pl

Arguments:

  \@fastq              The list of fastq files
  $expectedGenomeSize  in nucleotides
  $settings            hash

Returns

  $properties:  hash ref with keys.  E.g.,
    File
    avgReadLength
    totalBases
    minReadLength
    maxReadLength
    avgQuality
    numReads
    coverage

=cut

sub readMetrics{
  my($fastqs, $expectedGenomeSize, $settings) = @_;
  my %metrics;

  my @fastq = sort{$a cmp $b} @$fastqs;

  my $tempdir = tempdir(basename($0).".readMetrics.XXXXXX",TMPDIR=>1,CLEANUP=>1);

  for(my $i=0; $i<@fastq; $i++){
    my $seqtkComp = "$tempdir/".basename($fastq[$i]).".seqtkcomp";
    my $seqtkFqchk= "$tempdir/".basename($fastq[$i]).".seqtkfqchk";
    command("seqtk comp $fastq[$i] > $seqtkComp");
    command("seqtk fqchk $fastq[$i]> $seqtkFqchk");

    my(@lengths,@A,@C,@G,@T);
    open(my $fh, $seqtkComp) or die "ERROR: could not read seqtk comp file $seqtkComp: $!";
    while(<$fh>){
      # chomp; # <- not including the last field and so no need to remove last newline char
      my($readName, $length, $A, $C, $G, $T) = split /\t/;
      push(@lengths, $length);
      push(@A, $A);
      push(@C, $C);
      push(@G, $G);
      push(@T, $T);
    }
    close $fh;

    # parse Fqchk to get average qual
    my $avgQual = -1;
    open(my $fh2, $seqtkFqchk) or die "ERROR: could not read seqtk fqchk file $seqtkFqchk: $!";
    while(<$fh2>){
      next if(!/^ALL/);
      chomp;
      my(undef, $bases, $A,$C,$G,$T,$N, $avgQ, $errQ, $low, $high) = split /\t/;
      $avgQual = $avgQ;
      last;
    }
    close $fh2;

    my $totalBases = sum(@lengths);
    my $numReads   = scalar(@lengths);
    my $coverage   = -1;
    if($expectedGenomeSize){
      $coverage = $totalBases / $expectedGenomeSize;
    }
    $metrics{$fastq[$i]} = {
      avgReadLength  => ($totalBases/$numReads),
      totalBases     => $totalBases,
      minReadLength  => min(@lengths),
      maxReadLength  => max(@lengths),
      avgQuality     => $avgQual,
      numReads       => $numReads,
      coverage       => $coverage,
    };
  }
  return \%metrics;
}

1;

