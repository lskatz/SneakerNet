SYNOPSIS
========

Any executable in this directory will be run when SneakerNet is executed.  If 
permissions are not set as executable, then the script will not be run as a
plugin.  A good way to test a new script, therefore, is to introduce it to this 
folder but not make it executable.

Each executable will be run with standard options and ARGV=="run directory"

How to make a plugin
====================

1. Place it in the SneakerNet.plugins folder.
2. Test the script with a run directory.
  1. The first positional argument must be the MiSeq run directory
  2. The script must accept the following flags with the following example 
     values (or no values)
     1. `--help`
     2. `--numcpus 1`
     3. `--debug`
     4. `--tempdir /some/directory`
     5. `--help`
  3. Add any desired soft-coded variables into config.bak/settings and config/settings
  4. If any files are added to the run directory, please add them to 
     `runDirectory/SneakerNet/customdirectory`, and add any results for the
     resulting email to `runDirectory/SneakerNet/forEmail`. Any files under
     this directory will be emailed with the SneakerNet email.
3. If everything runs properly, chmod the script for the sequencing user, so
   that it is executable.