SYNOPSIS
========

Any executable in this directory will be run when SneakerNet is executed.  If 
permissions are not set as executable, then the script will not be run as a
plugin.  A good way to test a new script, therefore, is to introduce it to this 
folder but don't make it executable.

Each executable will be run with standard options and ARGV=="run directory"

How to make a plugin
====================

1. Place it in the SneakerNet.plugins folder.
2. Test the script with a run directory.
  1. The first positional argument must be the MiSeq run directory
  2. The script must accept the following flags with the following example 
     values
     1. `--help`
     2. `--numcpus 1`
     3. `--debug`
3. If everything runs properly, chmod the script for the sequencing user, so
   that it is executable.