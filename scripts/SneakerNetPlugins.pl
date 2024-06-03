#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use Cwd qw/getcwd realpath/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;
use Config::Simple;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;
my $cwd = getcwd();

# All logging will go to a file, which will end up as $run/SneakerNet.txt.
# The link to the log will be emailed in the email plugin.

exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(dry-run keep-going tempdir=s help version numcpus=i no|skip=s@ email! force! workflow=s)) or die $!;

  if($$settings{version}){
    print "SneakerNet $SneakerNet::VERSION\n";
    return 0;
  }

  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{email}//=1;
  $$settings{no}   //=[];
  if($$settings{tempdir}){
    mkdir $$settings{tempdir};
  }

  # Have a not-so-hidden Konami code Easter egg for fun
  if(join(" ", @ARGV) =~ /up up down down left right left right B A( start)?/i){
    logmsg "The great bioinformatics god in the sky blesses this run.";
    exit(0);
  }

  # The rest of the positional arguments are the run director(ies)
  my @dir=@ARGV;

  for my $dir(@dir){
    # Make the basic SneakerNet folders
    for("$dir/SneakerNet", "$dir/SneakerNet/forEmail"){
      mkdir $_;
    }

    # Determine the workflow
    my $cfg = new Config::Simple();
    my %vars;
    if($cfg->read("$dir/snok.txt")){
      %vars = $cfg->vars();
    }
    # Workflow is either defined in CLI settings, snok.txt, or is set to default.
    my $workflow = $$settings{workflow} || $vars{'default.workflow'} || "default";

    my $exe = $$settings{"plugins.$workflow"};
    if(!defined($exe) || ref($exe) ne 'ARRAY'){
      die "ERROR: workflow in snok.txt was defined as $workflow, but this workflow is not defined in config/plugins.conf";
    }

    # Get any additional plugins requested
    my $additionalPlugins = $vars{"default.additional_plugins"} || [];
    $additionalPlugins = (ref($additionalPlugins) eq 'ARRAY')?$additionalPlugins:[$additionalPlugins];
    # Add these plugins to the workflow
    push(@$exe, @$additionalPlugins);

    # Should we skip any plugins?
    # First, email since that is its own parameter
    if(!$$settings{email}){
      $exe = [grep {$_ !~ /email/} @$exe];
    }
    # Additional skips
    for my $skipPlugin(@{ $$settings{no} }){
      $exe = [grep {$_ !~ /$skipPlugin/i} @$exe];
    }
    
    # Run all plugins
    chdir($dir) or die "ERROR: could not change to directory $dir: $!";
    logmsg "Changing directory to $dir";
    for my $e(@$exe){

      my $command="$RealBin/../SneakerNet.plugins/$e . --numcpus $$settings{numcpus}";
      $command.=" --force" if($$settings{force});
      $command.=" --tempdir $$settings{tempdir}" if($$settings{tempdir});

      if($$settings{'dry-run'}){
        logmsg $command;
        next;
      }
      eval{
        command($command);
      };
      if($@){
        logmsg "ERROR with plugin $e: $@";
        if($$settings{'keep-going'}){
          logmsg "  ... however, --keep-going was specified and I will ignore that.";
          next;
        }
        die "  ... crashing because --keep-going was not specified";
      }
    }
    chdir($cwd) or die "ERROR: could not change back to directory $cwd: $!"; # go back to original directory
    logmsg "Changing directories from $dir, back to $cwd";
  }
  
  return 0;
}

sub usage{
  "$0: runs all SneakerNet plugins on a run directory
  Usage: $0 dir [dir2...]
  --noemail     Do not send an email at the end.
  --no      ''  Specify a regex of which plugins to skip.
                Can provide multiple instances.  E.g.,
                --no email --no transfer --no save
                Case insensitive.
  --dry-run     Just print the plugin commands that would have been run
  --keep-going  If a plugin has an error, move onto the next anyway
  --numcpus 1
  --tempdir ''  Force a temporary directory path to each plugin
  --force
  --version     Print SneakerNet version and exit
  --workflow    Which workflow under plugins.conf should we follow?
                If not specified, will look at snok.txt.
                If not snok.txt, will use 'default'
  "
}
