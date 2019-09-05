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

| Plugin                         | description |
|:-------------------------------|:------------|
|sn_SalmID.pl                    | Salmonella subspecies identification |
|sn_passfail.pl                  | Table of pass/fail for each sample   |
|sn_iontorrent_assembleAll.pl    | Assembly for ion torrent data        |
|addReadMetrics.pl               | Raw read metrics                     |
|sn_helloWorld.py                | Example plugin in Python             |
|baseBalance.pl                  | Dividing all As by Ts and all Cs by Gs to see if we get a ratio of 1 for each|
|sn_mlst.pl                      | Runs 7-gene MLST on assemblies       |
|transferFilesToRemoteComputers.pl|Transfers files to a remote computer |
|sn_detectContamination.pl       | Detects potential contamination by kmer counting|
|emailWhoever.pl                 | Emails all results                   |
|sn_helloWorld.pl                | Example plugin in Perl               |
|sn_detectContamination-mlst.pl  | Runs 7-gene MLST on raw reads, checking for abnormal number of alleles |
|sn_iontorrent_parseSampleSheet.pl|Turns the sample sheet for ion torrent into SneakerNet format |
|sn_immediateStatus.pl           | Emails an immediate report           |
|guessTaxon.pl                   | Runs metagenomics classifier to guess the taxon for a sample |
|assembleAll.pl                  | Assembles Illumina data              |
|sn_parseSampleSheet.pl          | Turns the sample sheet for Illumina into SneakerNet format |
|sn_report.pl                    | Creates an HTML report from all other plugins |

## Workflows

Workflows place the plugins in a particular order for a given run.

### Default

MiSeq, HiSeq, MiniSeq, ...

### Ion Torrent

