#!/usr/bin/env perl

use strict;
use warnings;
use File::Copy qw/cp mv/;
use File::Basename qw/basename/;
 
use 5.010;
 
use ExtUtils::MakeMaker;

# Copy the config.bak directory
mkdir "config";
for my $file(glob("config.bak/*")){
  my $target = "config/".basename($file);
  next if(-e $target);
  cp($file, $target) or die $!;
}

# Some help from
# http://blogs.perl.org/users/michal_wojciechowski/2011/11/github-friendly-readme-files-with-extutils-makemaker-and-module-build.html
my $preop =
  'perldoc -uT $(VERSION_FROM) | tee $(DISTVNAME)/README.pod > README.pod;' .
  'pod2text README.pod | tee $(DISTVNAME)/README > README;' .
  #'mkdir -pv config' .
  #'cp -v config.bak/ config/' .
  '';

WriteMakefile(
  NAME         => 'SneakerNet',
  VERSION_FROM => 'lib/perl5/SneakerNet.pm',
  #ABSTRACT_FROM => 'lib/perl5/SneakerNet.pm',
  AUTHOR => q{Lee S. Katz <gzu2@cdc.gov>},
  EXE_FILES    => [
    glob("scripts/*.pl"),
    glob("SneakerNet.plugins/*.pl"),
    glob("SneakerNet.plugins/*.py"),
    glob("SneakerNet.plugins/*.sh"),
  ],
  PREREQ_PM    => {
    # Core modules
    'File::Basename' => 0,
    'Data::Dumper'   => 0,
    'List::Util'     => 0,
    'Exporter'       => 0,
    'TAP::Parser::YAMLish::Writer' => '3.44',
    # Threads modules
    'threads'        => 0,
    'threads::shared'=> 0,
    'Thread::Queue'  => 0,
    # Not core (I think?)
    'Statistics::Descriptive' => 0,
    'Config::Simple' => '4.58',
    'File::Slurp'    => 0,
    #'Bio::Kmer'      => 0, # Bio::Kmer currently requires Jellyfish v2 but Kraken1 requires Jellyfish v1
    'Bio::SeqIO'     => 0,
    'Text::Levenshtein'    => 0.15,
    'Bio::FeatureIO' => 0,
    # MLST module
    'Moo'            => 0,
    'List::MoreUtils'=> 0,
    'JSON'           => 0,
    'GD'             => 0,
    'GD::Graph::lines'=>0,
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

