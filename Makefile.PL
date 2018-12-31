use strict;
use warnings;
 
use 5.010;
 
use ExtUtils::MakeMaker;

my $VERSION = `grep -i version config.bak/settings.conf  | cut -f 2`;
chomp $VERSION;

# Some help from
# http://blogs.perl.org/users/michal_wojciechowski/2011/11/github-friendly-readme-files-with-extutils-makemaker-and-module-build.html
my $preop =
  #'sed -i.bak "s/.5f/.10f/g" bin/quicktree-2.3/src/tree.c bin/quicktree-2.3/src/distancemat.c;' .
  #'make -C bin/quicktree-2.3;' .
  'perldoc -uT $(VERSION_FROM) | tee $(DISTVNAME)/README.pod > README.pod;' .
  'pod2text README.pod | tee $(DISTVNAME)/README > README';

WriteMakefile1(
  NAME         => 'SneakerNet',
  VERSION      => $VERSION,
  ABSTRACT_FROM => 'README.md',
  AUTHOR => q{Lee S. Katz <gzu2@cdc.gov>},
  EXE_FILES    => [
    "scripts/SneakerNet.ionTorrent.ro.pl",
    "scripts/SneakerNet.pl",
    "scripts/SneakerNetPlugins.pl",
    "scripts/SneakerNet.roRun.pl",
    "SneakerNet.plugins/addReadMetrics.pl",
    "SneakerNet.plugins/assembleAll.pl",
    "SneakerNet.plugins/baseBalance.pl",
    "SneakerNet.plugins/emailWhoever.pl",
    "SneakerNet.plugins/guessTaxon.pl",
    "SneakerNet.plugins/read_length_quality_stats_fastq.py",
    "SneakerNet.plugins/sn_detectContamination.pl",
    "SneakerNet.plugins/sn_mlst.pl",
    "SneakerNet.plugins/sn_parseSampleSheet.pl",
    "SneakerNet.plugins/sn_passfail.pl",
    "SneakerNet.plugins/sn_SalmID.pl",
    "SneakerNet.plugins/transferFilesToRemoteComputers.pl",
  ],
  PREREQ_PM    => {
    # Core modules
    'File::Basename' => 0,
    'Data::Dumper'   => 0,
    'List::Util'     => 0,
    'Exporter'       => 0,
    # Threads modules
    'threads'        => 0,
    'threads::shared'=> 0,
    'Thread::Queue'  => 0,
    # Not core (I think?)
    'Statistics::Descriptive' => 0
    'Config::Simple' => 0,
    'File::Slurp'    => 0,
    'Bio::Kmer'      => 0,
    'Text::Fuzzy'    => 0,
    'Email::Stuffer' => 0,
    'Bio::SeqIO'     => 0,
    'Bio::FeatureIO' => 0,
    # MLST module
    'Moo'            => 0,
    'List::MoreUtils'=> 0,
    'JSON'           => 0,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
  MIN_PERL_VERSION => '5.10.0',
  LICENSE          => "gpl_3",
  dist=>{
    COMPRESS => "gzip -9f",
    SUFFIX   => "gz",
    PREOP    => $preop,
  },
  META_MERGE => {
    'meta-spec' => { version => 2 },
    resources => {
      repository => {
        type => 'git',
        url => 'git://github.com/lskatz/SneakerNet.git',
        web => 'https://github.com/lskatz/SneakerNet',
      },
    },
  },
);

sub WriteMakefile1 {  #Compatibility code for old versions of EU::MM. Written by Alexandr Ciornii, version 2. Added by eumm-upgrade.
  my %params=@_;
  my $eumm_version=$ExtUtils::MakeMaker::VERSION;
  $eumm_version=eval $eumm_version;
  die "EXTRA_META is deprecated" if exists $params{EXTRA_META};
  die "License not specified" if not exists $params{LICENSE};
  if ($params{AUTHOR} and ref($params{AUTHOR}) eq 'ARRAY' and $eumm_version < 6.5705) {
    $params{META_ADD}->{author}=$params{AUTHOR};
    $params{AUTHOR}=join(', ',@{$params{AUTHOR}});
  }
  if ($params{TEST_REQUIRES} and $eumm_version < 6.64) {
    $params{BUILD_REQUIRES}={ %{$params{BUILD_REQUIRES} || {}} , %{$params{TEST_REQUIRES}} };
    delete $params{TEST_REQUIRES};
  }
  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
    #EUMM 6.5502 has problems with BUILD_REQUIRES
    $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
    delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;

  WriteMakefile(%params);
}