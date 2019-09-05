# SYNOPSIS

SneakerNet is a set of plugins. Each plugin runs a distinct analysis
on a MiSeq run and accepts specific flags when it is run.
Therefore, any given plugin can run independently of the others (aside
from any prerequisite files, e.g., genome assemblies from a previous
plugin).

## Command line

Each plugin can accept the following options. The first positional
parameter must be the SneakerNet run.

|Flag|Default value|description|
|:---|:------------|:-----------|
|`--help`|         |generate a help menu|
|`--numcpus`|     1|Parallelization|
|`--debug`|        |generate more messages or any other debugging|
|`--tempdir`|automatically generated, e.g., with `File::Temp` or `mktemp`|Where temporary files are located|
|`--force`|        |This is loosely defined but can be used for many things like overwriting output files|
|`--version`|      |Print a version in the format of X.Y or X.Y.Z|
|`--citation`|     | Print a citation statement. | 


## Catalog

Except for the legacy plugins, all plugins are prefixed with `sn_`.
The plugins are not specific to any one language, although the majority
are in Perl.

sn_SalmID.pl
sn_passfail.pl
sn_iontorrent_assembleAll.pl
addReadMetrics.pl
sn_helloWorld.py
baseBalance.pl
sn_mlst.pl
transferFilesToRemoteComputers.pl
sn_detectContamination.pl
emailWhoever.pl
sn_helloWorld.pl
sn_detectContamination-mlst.pl
sn_iontorrent_parseSampleSheet.pl
sn_immediateStatus.pl
guessTaxon.pl
assembleAll.pl
sn_parseSampleSheet.pl
sn_report.pl

## Workflows

Workflows place the plugins in a particular order for a given run.

### Default

MiSeq, HiSeq, MiniSeq, ...

### Ion Torrent

