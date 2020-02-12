# SYNOPSIS

SneakerNet is a set of plugins. Each plugin runs a distinct analysis
on a MiSeq run and accepts specific flags when it is run.
Therefore, any given plugin can run independently of the others (aside
from any prerequisite files, e.g., genome assemblies from a previous
plugin).

# Workflows

SneakerNet workflows define a particular order for the plugins to run.
They help resolve dependencies like ensuring that genome assemblies
are present before analyzed or enforcing that a report is generated only
after all plugins have created their outputs.

## Default

This workflow runs most plugins and assumes that you have some flavor
of Illumina (MiSeq, HiSeq, MiniSeq).

## Ion Torrent

This workflow runs plugins designed for ion torrent.

## Metagenomics

For metagenomics runs. 

# Command line

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
|`--check-dependencies`|     | check all executable dependencies | 


# Catalog

Except for the legacy plugins, all plugins are prefixed with `sn_`.
The plugins are not specific to any one language, although the majority
are in Perl.

Contributions are welcome for the following plugin documents.
Some links might not be valid at this time.

| Plugin                         | description |
|:-------------------------------|:------------|
|[sn_SalmID.pl](plugins/sn_SalmID.pl.md)          | Salmonella subspecies identification |
|[sn_passfail.pl](plugins/sn_passfail.pl.md)      | Table of pass/fail for each sample   |
|[sn_iontorrent_assembleAll.pl](plugins/sn_iontorrent_assembleAll.pl.md)    | Assembly for ion torrent data        |
|[addReadMetrics.pl](plugins/addReadMetrics.pl.md)| Raw read metrics                     |
|[sn_helloWorld.pl](plugins/sn_helloWorld.pl.md)               | Example plugin in Perl               |
|[sn_helloWorld.sh](plugins/sn_helloWorld.sh.md)               | Example plugin in Bash               |
|[sn_helloWorld.py](plugins/sn_helloWorld.py.md)  | Example plugin in Python             |
|[baseBalance.pl](plugins/baseBalance.pl.md)      | Dividing all As by Ts and all Cs by Gs to see if we get a ratio of 1 for each|
|[sn_mlst.pl](plugins/sn_mlst.pl.md)              | Runs 7-gene MLST on assemblies       |
|[transferFilesToRemoteComputers.pl](plugins/transferFilesToRemoteComputers.pl.md)|Transfers files to a remote computer |
|[sn_detectContamination.pl](plugins/sn_detectContamination.pl.md)       | Detects potential contamination by kmer counting|
|[emailWhoever.pl](plugins/emailWhoever.pl.md)                 | Emails all results                   |
|[sn_detectContamination-mlst.pl](plugins/sn_detectContamination-mlst.pl.md)  | Runs 7-gene MLST on raw reads, checking for abnormal number of alleles |
|[sn_iontorrent_parseSampleSheet.pl](plugins/sn_iontorrent_parseSampleSheet.pl.md)|Turns the sample sheet for ion torrent into SneakerNet format |
|[sn_immediateStatus.pl](plugins/sn_immediateStatus.pl.md)           | Emails an immediate report           |
|[guessTaxon.pl](plugins/guessTaxon.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample |
|[sn_detectContamination-kraken.pl](plugins/sn_detectContamination-kraken.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample and list at most a single major contaminant |
|[assembleAll.pl](plugins/assembleAll.pl.md)                  | Assembles Illumina data              |
|[sn_crypto_assembleAll.pl](plugins/sn_crypto_assembleAll.pl.md)                  | Assembles Illumina data for Cryptosporidium        |
|[sn_parseSampleSheet.pl](plugins/sn_parseSampleSheet.pl.md)          | Turns the sample sheet for Illumina into SneakerNet format |
|[sn_report.pl](plugins/sn_report.pl.md)                    | Creates an HTML report from all other plugins |
|[sn_saveFailedGenomes.pl](plugins/sn_saveFailedGenomes.pl.md)                    | Saves genomes into the destination folder, into a QC_Fails subfolder|

