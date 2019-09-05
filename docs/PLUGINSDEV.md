# SYNOPSIS

These are instructions on how to make a SneakerNet plugin.  It does not
matter which language the plugin is coded in.  All that matters is that
it is executable by the SneakerNet user and that it can accept certain
parameters.

# EXAMPLES

Two plugin scripts have been provided in the plugins directory.
`sn_helloWorld.pl` is coded in Perl, and 
`sn_helloWorld.py` is coded in Python.

# How to make a plugin

Two major steps are described below for making a plugin.

## Create the script.
1. The first positional argument must be the MiSeq run directory
2. The script must accept the following flags with the following example 
   values (or no values). The flags must adhere to original meaning behind
   the flags. For example, `--help` must provide a help menu.
   The script does not necessarily need to _use_ these flags however.
3. Add any desired soft-coded variables such as the location of a blast database
   into `config.bak/settings` and `config/settings`.
   `config.bak` is tracked with `git`, and `config` is created upon
   SneakerNet installation.
4. Inputs: the plugin can and should read sample names and properties
   from `samples.tsv`.
5. If the plugin generates any files, please organize them into 
   `runDirectory/SneakerNet/customdirectory` (where `customdirectory`
   is a name of your choice), and add any results for the
   resulting email to `runDirectory/SneakerNet/forEmail`. Any files under
   this directory will be emailed with the SneakerNet email.
6. List your results in `runDirectory/SneakerNet/properties.tsv`. This will 
   ensure that your results will be in the email report.
   1. `properties.tsv` is three columns: plugin name, key, value
   2. You should include `version`
     
|Flag|Default value|description|
|:---|:------------|:-----------|
|`--help`|         |generate a help menu|
|`--numcpus`|     1|Parallelization|
|`--debug`|        |generate more messages or any other debugging|
|`--tempdir`|automatically generated, e.g., with `File::Temp` or `mktemp`|Where temporary files are located|
|`--force`|        |This is loosely defined but can be used for many things like overwriting output files|
|`--version`|      |Print a version in the format of X.Y or X.Y.Z|
|`--citation`|     | Print a citation statement to give yourself credit. | 

## Activate the script as a plugin

1. Place it in the SneakerNet.plugins folder
2. chmod the script to be executable
3. Add the plugin to the list of plugins in `config.bak/plugins` and `config/plugins` 

## Unit tests

Create a file under `t/` and copy a previous example unit test from that directory.

To run all unit tests:

    perl Makefile.PL
    make
    make test
