# SYNOPSIS

SneakerNet is a set of plugins. Each plugin runs a distinct analysis
on a MiSeq run and accepts specific flags when it is run.
Therefore, any given plugin can run independently of the others (aside
from any prerequisite files, e.g., genome assemblies from a previous
plugin).

Quick table of contents:
- [Workflows](#workflows)
- [Which plugins are running?](#workflows-and-their-plugins)
- [Plugins on the command line](#command-line)
- [Plugins catalog](#catalog)

# Workflows

SneakerNet workflows define a particular order for the plugins to run.
They help resolve dependencies like ensuring that genome assemblies
are present before analyzed or enforcing that a report is generated only
after all plugins have created their outputs.

Workflows are defined in [plugins.conf](../config.bak/plugins.conf).

The exact order of plugins for all workflows can be found by running the command `SneakerNet.checkdeps.pl --list`.

To make your own custuom workflow, edit the file under `config/plugins.conf`.
Plugins are run in the order specified for any given workflow.  For example:

    default = pluginA.pl, pluginB.pl, pluginZ.pl

In this example, in the default workflow, 
if `pluginB` has a dependency on `pluginZ`, you might want to change the order
so that `pluginB` runs last.

    default = pluginA.pl, pluginZ.pl, pluginB.pl

## Default

This workflow runs most plugins and assumes that you have some flavor
of Illumina (MiSeq, HiSeq, MiniSeq).

## Ion Torrent

This workflow runs plugins designed for ion torrent.

## Metagenomics

For metagenomics runs. 

## Assembly

For assembly-only runs (ie, only assemblies and not raw reads in the folder).

## sarscov2

For running the SARS-CoV-2 workflow. Plugin(s) are prefixed with `sn_sars_`.

# Workflows and their plugins

Workflows in SneakerNet are defined by which plugins are run and in which order.
You can run `SneakerNet.checkdeps.pl --list` to see which plugins are run and in which order, for any workflow.
This script pulls from `config/plugins.conf`. 
Below is the output for SneakerNet version 0.11.2.
By default, only `default` will be run in a SneakerNet analysis, but other workflows are available.
The pseudo-workflow `all` is an alphabetical listing of all available plugins.

    $ SneakerNet.checkdeps.pl --list
        all
        addReadMetrics.pl, assembleAll.pl, baseBalance.pl, emailWhoever.pl, sn_assemblyWorkflow_init.pl, sn_crypto_assembleAll.pl, sn_crypto_gp60.pl, sn_detectContamination-kraken.pl, sn_detectContamination-mlst.pl, sn_detectContamination.pl, sn_helloWorld.pl, sn_helloWorld.py, sn_helloWorld.sh, sn_immediateStatus.pl, sn_iontorrent_assembleAll.pl, sn_kraken-metagenomics.pl, sn_kraken.pl, sn_mlst.pl, sn_mlst-wg.pl, sn_parseSampleSheet.pl, sn_passfail.pl, sn_report.pl, sn_SalmID.pl, sn_saveFailedGenomes.pl, sn_staramr.pl, transferFilesToRemoteComputers.pl

        assembly
        sn_assemblyWorkflow_init.pl, sn_mlst.pl, sn_staramr.pl, sn_passfail.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_report.pl, emailWhoever.pl

        cryptosporidium
        sn_parseSampleSheet.pl, addReadMetrics.pl, sn_crypto_assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_passfail.pl, transferFilesToRemoteComputers.pl, emailWhoever.pl

        default
        sn_parseSampleSheet.pl, addReadMetrics.pl, assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_detectContamination-mlst.pl, baseBalance.pl, sn_staramr.pl, sn_passfail.pl, transferFilesToRemoteComputers.pl, sn_report.pl, emailWhoever.pl

        iontorrent
        addReadMetrics.pl, sn_iontorrent_assembleAll.pl, sn_mlst.pl, sn_kraken.pl, sn_detectContamination-kraken.pl, sn_passfail.pl, sn_staramr.pl, transferFilesToRemoteComputers.pl, emailWhoever.pl

        metagenomics
        sn_parseSampleSheet.pl, addReadMetrics.pl, sn_kraken.pl, sn_kraken-metagenomics.pl, sn_passfail.pl, sn_report.pl

![Default workflow](/docs/images/defaultworkflow.png)

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
|`--check-dependencies`|     | check all executable dependencies. Print executable dependencies to stdout and version information to stderr. Run `SneakerNet.checkdeps.pl` to check dependencies on all plugins. | 


# Catalog

Except for the legacy plugins, all plugins are prefixed with `sn_`.
The plugins are not specific to any one language, although the majority
are in Perl.

Contributions are welcome for the following plugin documents.

| Plugin                         | description |
|:-------------------------------|:------------|
|[sn_SalmID.pl](plugins/sn_SalmID.pl.md)          | Salmonella subspecies identification |
|[sn_staramr.pl](plugins/sn_staramr.pl.md)          | staramr antimicrobial resistance determinant analysis |
|[sn_passfail.pl](plugins/sn_passfail.pl.md)      | Table of pass/fail for each sample   |
|[sn_iontorrent_assembleAll.pl](plugins/sn_iontorrent_assembleAll.pl.md)    | Assembly for ion torrent data        |
|[addReadMetrics.pl](plugins/addReadMetrics.pl.md)| Raw read metrics                     |
|[sn_helloWorld.pl](plugins/sn_helloWorld.pl.md)               | Example plugin in Perl               |
|[sn_helloWorld.sh](plugins/sn_helloWorld.sh.md)               | Example plugin in Bash               |
|[sn_helloWorld.py](plugins/sn_helloWorld.py.md)  | Example plugin in Python             |
|[baseBalance.pl](plugins/baseBalance.pl.md)      | Dividing all As by Ts and all Cs by Gs to see if we get a ratio of 1 for each|
|[sn_mlst.pl](plugins/sn_mlst.pl.md)              | Runs 7-gene MLST on assemblies       |
|[sn_mlst-wg.pl](plugins/sn_mlst-wg.pl.md)              | Runs whole-genome MLST on assemblies       |
|[transferFilesToRemoteComputers.pl](plugins/transferFilesToRemoteComputers.pl.md)|Transfers files to a remote computer |
|[sn_detectContamination.pl](plugins/sn_detectContamination.pl.md)       | Detects potential contamination by kmer counting|
|[emailWhoever.pl](plugins/emailWhoever.pl.md)                 | Emails all results                   |
|[sn_detectContamination-mlst.pl](plugins/sn_detectContamination-mlst.pl.md)  | Runs 7-gene MLST on raw reads, checking for abnormal number of alleles |
|[sn_iontorrent_parseSampleSheet.pl](plugins/sn_iontorrent_parseSampleSheet.pl.md)|Turns the sample sheet for ion torrent into SneakerNet format |
|[sn_immediateStatus.pl](plugins/sn_immediateStatus.pl.md)           | Emails an immediate report           |
|[guessTaxon.pl](plugins/guessTaxon.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample |
|[sn_kraken.pl](plugins/sn_kraken.pl.md)                   | Runs metagenomics classifier on raw reads, or on assemblies if reads are not present. No secondary analysis is performed by this exact plugin. E.g., `sn_detectContamination-kraken.pl`. |
|[sn_detectContamination-kraken.pl](plugins/sn_detectContamination-kraken.pl.md)                   | Runs metagenomics classifier to guess the taxon for a sample and list at most a single major contaminant |
|[sn_kraken-metagenomics.pl](plugins/sn_kraken-metagenomics.pl.md)                   | Analyzes kraken results for a metagenomics sample |
|[assembleAll.pl](plugins/assembleAll.pl.md)                  | Assembles Illumina data              |
|[sn_assemblyWorkflow_init.pl](plugins/sn_assemblyWorkflow_init.pl.md)                  | For workflows that only have assembly data. Initializes the workflow so that other plugins can function properly. |
|[sn_crypto_assembleAll.pl](plugins/sn_crypto_assembleAll.pl.md)                  | Assembles Illumina data for Cryptosporidium        |
|[sn_crypto_gp60.pl](plugins/sn_crypto_gp60.pl.md)                  | Provides the gp60 profile for Cryptosporidium        |
|[sn_parseSampleSheet.pl](plugins/sn_parseSampleSheet.pl.md)          | Turns the sample sheet for Illumina into SneakerNet format |
|[sn_report.pl](plugins/sn_report.pl.md)                    | Creates an HTML report from all other plugins |
|[sn_sarscov2_assembleAll.pl](plugins/sn_sarscov2_assembleAll.pl.md)              | Runs assembly for SARS-CoV-2 amplicon-based genomes |
|[sn_assembleAll_reference.pl.md](plugins/sn_assembleAll_reference.pl.md)              | Runs reference assembly |
|[sn_saveFailedGenomes.pl](plugins/sn_saveFailedGenomes.pl.md)                    | Saves genomes into the destination folder, into a QC_Fails subfolder|

