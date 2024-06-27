# SneakerNet Output

SneakerNet creates a subdirectory `SneakerNet/` in the input folder.
In this document, this subfolder will be referred to as `SneakerNet/`.
This document describes those files.

## Files under `SneakerNet/`

### properties.txt

This file records special values from each plugin in the format of
plugin name, key, value.
Duplicate plugin name and key combinations are allowed; however,
the last plugin name and key combination overrides the earlier instances.
Therefore for example, a warning from a plugin will be "erased"
if the plugin is run again and it records a value describing zero warnings.
There are some special keys that are treated differently:

* table - the relative path to a `.tsv` file from the main folder that describes a summary table from a plugin, e.g., `SneakerNet/forEmail/kraken.tsv`. This can also be a csv file.
* version - the version number of the plugin that was run in the format of simple integers separated by periods, e.g., `2.3.1`
* XX-version - the version number of some software used in this specific plugin
* warnings - a summary of warnings from the plugin, e.g., `5 samples do not have fastq files`. Use "0" for no warnings.
* errors - a summary of errors from the plugin, e.g., `Error reading summary table`. Use "0" for no errors.
* date - date of analysis in YYYY-MM-DD format
* time - time of analysis in HH:MM:SS format
* other - other keys are allowed but not not specifically described here, e.g., `resfinder_gene_drug_version` for the staramr plugin. For other special keys, please see the individual plugin documentation pages.

### forEmail/passfail.tsv

This file records sample names and categories on which
a given sample can pass or fail. Each value for a category
is a simple integer.

* 1 - This sample failed on this category
* 0 - This sample passed on this category
* -1 - It is unknown whether this sample passed. For example, if the taxon is unknown and therefore the coverage is unknown.

### Other files under `forEmail/`

Please see individual plugin documentation pages on other files

## Subfolders

There are some special subfolders in `SneakerNet/`.

### forEmail

`SneakerNet/forEmail/` is the directory that all email attatchments live.
Plugins will create summary files, usually tables, and then they will get attached to the report email.

### assemblies

`SneakerNet/assemblies` is the directory that all genome assemblies and gene predictions live.
For each sample, e.g., `sample1`, there is a subfolder, e.g., `SneakerNet/assemblies/sample1`.
There is at least one `.fasta` and `.gbk` file in each sample subfolder.
Each file is named after the basic assembly method, e.g., `SneakerNet/assemblies/sample1/sample1.skesa.fasta`.

### kraken

`SneakerNet/kraken` contains a subfolder for each sample name,
e.g., `SneakerNet/kraken/sample1`.
In each sample folder, there are
kraken output files, or output files compatible with Kraken output files.
These files are:

* kraken.report - a tab-separated file whose values have optional space padding. Its fields are
  * percent of reads
  * number of reads
  * number of reads specific to this exact taxon
  * rank - one letter representation, e.g., `S` for species
  * taxid
  * Name for the taxon
* kraken.taxonomy - a tab-separated file whose values are
  * number of reads specific to this exact taxon
  * taxonomy lineage separated by tabs, starting with `root` and ending with the genus/species name
* report.html - the Krona html file for visualization

`kraken.out` is _not_ always present due to space limitations but
would describe the taxonomic classification for each read.

### Other

Other plugins can create subfolders. For example, the plugin `baseBalance.pl`
creates a subfolder `SneakerNet/baseBalance`.
Please see individual plugin documentation pages for more information.


