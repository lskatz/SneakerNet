#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;

use Test::More tests => 2;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use SneakerNet;

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
      diag "Skipping plugin $file - not executable";
      next;
    }
    # TODO add long options to the shell script and then remove this if block
    if($path =~ /\.sh$/){
      diag "Skipping $file - shell script";
      next;
    }
    push(@plugin, $path);
  }
}
closedir($dh);

subtest 'All plugins listed in PLUGINS.md' => sub{
  plan tests => scalar(@plugin);

  my $plugins_md = "$RealBin/../docs/PLUGINS.md";
  diag "Reading from $plugins_md";
  open(my $fh, $plugins_md) or BAIL_OUT("ERROR: could not read $plugins_md: $!");
  my @plugins_md_content = <$fh>;
  close $fh;
  #note "@plugins_md_content";

  for my $plugin(@plugin){
    my $exe = basename($plugin);
    my @documentationLine = grep {/$exe/i} @plugins_md_content;
    cmp_ok(scalar(@documentationLine), '>', 0, "$exe mentioned at least once in $plugins_md");
  }
};

subtest 'All plugins have a documentation page' => sub{
  plan tests => 1;

  pass("This test only shows warnings");

  for my $plugin(@plugin){
    my $exe = basename($plugin);
    my $md = "$RealBin/../docs/plugins/$exe.md";
    if(! -e $md){
      note "WARNING: could not find a documentation page for $exe under $md";
    }
  }
};

