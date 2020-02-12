#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet qw/logmsg/;

$ENV{PATH}="$RealBin/../scripts:$RealBin/../SneakerNet.plugins:$ENV{PATH}";

# Gather all plugins
my @plugin;
my $pluginsDir = "$RealBin/../SneakerNet.plugins";
opendir(my $dh, $pluginsDir) or die "ERROR reading directory $pluginsDir: $!";
while(my $file = readdir($dh)){
  next if($file=~/\.tmp$/); # Don't even worry about the tmp files that might be present
  my $path = "$pluginsDir/$file";
  if(-f $path){
    if(! -x $path){
      logmsg "Skipping plugin $file - not executable";
      next;
    }
    # TODO add long options to the shell script and then remove this if block
    if($path =~ /\.sh$/){
      logmsg "Skipping $file - shell script";
      next;
    }
    push(@plugin, $path);
  }
}
closedir($dh);

my $err_count = 0;
for my $path(@plugin){
  system("$path --check-dependencies");
  my $is_error = !! $?;
  $err_count += $is_error;
  #die if $?;
}

logmsg "Number of plugins with missing dependencies: $err_count";

