#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use Getopt::Long qw/GetOptions/;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/logmsg readConfig/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";

our $VERSION = 1;

exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings, qw(list help)) or die $!;
  $$settings{help} && usage();

  if($$settings{list}){
    # List all workflows and their plugins, and then exit
    my $workflows = $$settings{obj}{'plugins.conf'}->param(-block=>'plugins');
    for my $w (sort keys %$workflows){
      print "$w\n";
      print join(", ", @{$$workflows{$w} })."\n\n";
    }
    return 0;
  }

  # Which workflows does the user want to talk about
  my @workflow = @ARGV;
  if(!@workflow){
    logmsg "No workflow given. Will use 'default'";
    $workflow[0] = 'default';
  }

  # Get all dependencies per workflow
  my %dep;
  for my $w(@workflow){
    my $dep = checkWorkflowDependencies($w, $settings);
    %dep = (%dep, %$dep);
  }

  # Print dependencies
  print "DEPENDENCIES THAT WERE CHECKED:\n";
  for my $d(sort keys %dep){
    print "$d\n";
  }

  return 0;
}

sub checkWorkflowDependencies{
  my ($workflow, $settings) = @_;

  my $pluginArr = $$settings{"plugins.$workflow"};
  my $pluginsDir = "$RealBin/../SneakerNet.plugins";

  if(ref($pluginArr) ne 'ARRAY'){
    die "ERROR: workflow $workflow is not defined";
  }
  if(@$pluginArr < 1){
    die "ERROR: workflow $workflow was defined but there are no plugins associated with it.";
  }

  my %dep;
  for my $plugin(@$pluginArr){
    my $path = "$pluginsDir/$plugin";
    my @dep = `$path --check-dependencies`;
    chomp(@dep);
    for my $d(@dep){
      $dep{$d} = 1;
    }
  }

  # A nice newline before stdout takes over
  print STDERR "\n";

  return \%dep;
}

sub usage{
  print "Checks dependencies of all plugins
  Usage: $0 [options] [workflow workflow2...]

  workflow    If given, checks only the plugins listed for
              a given workflow as defined in plugins.conf.
              Multiple workflows can be given.
              If not given, the default workflow will be used.
  --list      Lists possible workflows
  --versions  Lists versions (not yet implemented)
  ";
  exit 0;
}

