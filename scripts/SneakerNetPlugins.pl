#!/usr/bin/env perl
# Moves folders from the dropbox-inbox pertaining to reads, and
# Runs read metrics on the directory

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/fileparse basename dirname/;

use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use SneakerNet qw/readConfig command logmsg/;

$ENV{PATH}="$ENV{PATH}:/opt/cg_pipeline/scripts";

local $0=fileparse $0;

# All logging will go to a file, which will end up as $run/SneakerNet.txt.
# The link to the log will be emailed in the email plugin.

exit(main());

sub main{
  my $settings=readConfig();
  GetOptions($settings,qw(help numcpus=i email! force!)) or die $!;
  die usage() if($$settings{help} || !@ARGV);
  $$settings{numcpus}||=1;
  $$settings{email}//=1;

  my @dir=@ARGV;

  for my $dir(@dir){
    # Make the basic SneakerNet folders
    for("$dir/SneakerNet", "$dir/SneakerNet/forEmail"){
      mkdir $_;
    }
    
    # ensure that the samplesheet can be parsed. This is
    # a prerequisite for all plugins.
    command("$FindBin::RealBin/../SneakerNet.plugins/sn_parseSampleSheet.pl $dir");

    # Run all plugins
    my @exe=@{ $$settings{'plugins.default'} };
    for my $exe(@exe){
      if(!$$settings{email} && $exe=~/emailWhoever.pl/){
        next;
      }

      my $command="$FindBin::RealBin/../SneakerNet.plugins/$exe $dir --numcpus $$settings{numcpus}";
      $command.=" --force" if($$settings{force});
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
  "
}
