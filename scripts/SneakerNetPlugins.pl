#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;
use Cwd qw/getcwd/;

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
  GetOptions($settings,qw(help numcpus=i email! force! workflow=s)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{email}//=1;

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

    if(!$$settings{email}){
      $exe = [grep {$_ !~ /email/} @$exe];
    }
    
    # Run all plugins
    chdir($dir) or die "ERROR: could not change to directory $dir: $!";
    for my $e(@$exe){
      my $command="$RealBin/../SneakerNet.plugins/$e . --numcpus $$settings{numcpus}";
      $command.=" --force" if($$settings{force});
      #print "$command\n\n";next;
      command($command);
    }
    chdir($cwd); # go back to original directory
  }
  
  return 0;
}

sub usage{
  "$0: runs all SneakerNet plugins on a run directory
  Usage: $0 dir [dir2...]
  --noemail     # Do not send an email at the end.
  --numcpus 1
  --force
  --workflow    Which workflow under plugins.conf should we follow?
                If not specified, will look at snok.txt.
                If not snok.txt, will use 'default'
  "
}
