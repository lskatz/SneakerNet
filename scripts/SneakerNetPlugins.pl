#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;

use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;
use Config::Simple;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

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
    my $workflow = "default";
    my $cfg = new Config::Simple();
    if($cfg->read("$dir/snok.txt")){
      my %vars = $cfg->vars();
      $workflow = $vars{'default.workflow'};
    }
    # But if the user specified the workflow, overwrite it
    if($$settings{workflow}){
      $workflow = $$settings{workflow};
    }

    my $exe = $$settings{"plugins.$workflow"};
    if(!defined($exe) || ref($exe) ne 'ARRAY'){
      die "ERROR: workflow in snok.txt was defined as $workflow, but this workflow is not defined in config/plugins.conf";
    }

    if(!$$settings{email}){
      $exe = [grep {$_ !~ /email/} @$exe];
    }
    
    # Run all plugins
    for my $e(@$exe){
      my $command="$RealBin/../SneakerNet.plugins/$e $dir --numcpus $$settings{numcpus}";
      $command.=" --force" if($$settings{force});
      #print "$command\n\n";next;
      command($command);
    }
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
