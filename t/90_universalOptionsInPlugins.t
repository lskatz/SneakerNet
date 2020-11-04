#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use List::Util qw/uniq/;

use Test::More tests => 4;

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

# sort
@plugin = sort @plugin;

# Travis-CI makes it super hard to install certain things so...
if($ENV{TRAVIS} || $ENV{GITHUB_ACTIONS}){
  # Remove plugins with kraken, taxonomy, assembly, mlst, ...
  # because it is difficult to get travis to download and
  # install it all.
  my $regex = qx/kraken|taxon|assembl|mlst|crypto|staramr|salm/;
  my @removed = grep { /$regex/i} @plugin;
  @plugin     = grep {!/$regex/i} @plugin;
  @plugin = sort @plugin;
  note "Removed some plugins b/c CI environment:";
  note "  => ".join(" ", @removed);
}
note "Testing these plugins: ". join(", ", @plugin);

# Do we have all plugins here?  At least one!
cmp_ok(scalar(@plugin), '>', 1, "Gathering all plugins. ".scalar(@plugin)." found.");

# Test whether --flagopt works for each script/flag combination
my @depsFail = ();
my @pluginsDepsFail = ();
subtest 'flagopt' => sub{
  for my $path(@plugin){
    my $file = basename($path);
    for my $flagOpt(qw(debug force version citation check-dependencies)){
      my $stdout = "";
      my $exit_code = eval{
        $stdout = `$path --$flagOpt --help 2>$path.tmp`;
        my $exit_code = $? >> 8;
        return $exit_code;
      };
      is($exit_code, 0, "Plugin $file given --$flagOpt");
      if($exit_code){
        my $newlineChar=$/;
        push(@depsFail, split(/$newlineChar/, $stdout));
        push(@pluginsDepsFail, $file);

        local $/ = undef;
        open(my $fh, "$path.tmp");
        my $content = <$fh>;
        close $fh;
        diag("Plugin $file failed with exit code $exit_code. REASON:\n$content\n");
      }
      unlink("$path.tmp");
    }
  }
};
if(@depsFail){
  diag "Dependencies from plugins that failed were: ".join(", ",sort(uniq(@depsFail)));
  diag "This does not mean that every single executable failed; it just means that when a dependency failed, it lists all the dependency of that plugin.";
  diag "Plugins that failed dependency checks were: ".join(", ", sort(@pluginsDepsFail));
}

# Test whether --flagopt str works for script/flag combinations
subtest 'flag str' => sub{
  for my $path(@plugin){
    my $file = basename($path);
    for my $flag(qw(tempdir)){

      # Test some strings but not too wonky
      for my $str(qw(foo bar)){
        my $stdout = "";
        my $exit_code = eval{
          $stdout = `$path --$flag $str --help 2>$path.tmp`;
          my $exit_code = $? >> 8;
          return $exit_code;
        };
        is($exit_code, 0, "Plugin $file given --$flag $str");
        if($exit_code){
          local $/ = undef;
          open(my $fh, "$path.tmp");
          my $content = <$fh>;
          close $fh;
          diag("Plugin $file failed with exit code $exit_code. REASON:\n$content\n");
        }
        rmdir($str); # be sure that the output of --tempdir doesn't make things difficult
        unlink("$path.tmp");
      }
    }
  }
};


# Test whether --flagopt int works for script/flag combinations
subtest 'flag int' => sub{
  for my $path(@plugin){
    my $file = basename($path);
    for my $flag(qw(numcpus)){

      # Test some wonky integers
      for my $int(0, 1, 13, 73, 144){
        my $stdout = "";
        my $exit_code = eval{
          $stdout = `$path --$flag $int --help 2>$path.tmp`;
          my $exit_code = $? >> 8;
          return $exit_code;
        };
        is($exit_code, 0, "Plugin $file given --$flag $int");
        if($exit_code){
          local $/ = undef;
          open(my $fh, "$path.tmp");
          my $content = <$fh>;
          close $fh;
          diag("Plugin $file failed with exit code $exit_code. REASON:\n$content\n");
        }
        unlink("$path.tmp");
      }
    }
  }
};

      
