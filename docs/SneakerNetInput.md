# SneakerNet Input folder

## Required

### \*.fastq.gz

Raw read files, both R1 and R2 if available. Single end raw reads are okay too.

### samples.tsv

This file can be created from SampleSheet.csv by the plugin [sn_parseSampleSheet.pl](../docs/plugins/sn_parseSampleSheet.pl.md), from `SneakerNet.ro.pl --createsamplesheet`, or manually.
It has three columns, separated by tab: 

* sample - Name for a sample
* info - semicolon-delimited information. Keys are separated by an equals sign. For example: `taxon=Vibrio;route=calcEngine` means that the taxon is Vibrio and the ultimate destination for the files are in calcEngine. Vibrio is the key for the taxon in `config/taxonProperties.conf` and will not have meaning if there is no entry.
* filenames - semicolon-delimited fastq files. For example, `./2018AW-0585_S318_L001_R1_001.fastq.gz;./2018AW-0585_S318_L001_R2_001.fastq.gz`

## Recommended

### snok.txt

This file stands for "SneakerNet Okay" which is an indicator that the run is ready to be analyzed.
The file is only necessary if a custom script needs to know whether the run is ready to be analyzed in an automated fashion.
At minimum, `snok.txt` is a zero-byte file.

However, this file can also have content in a key/value format separated by an equals sign.
The values can be comma-separated. Whitespace around each value is ignored.  For example:

    emails   = nobody@gatech.edu, nobody@cdc.gov, noreply@example.com
    workflow = default

The current keys are `emails` and `workflow`.
Emails describe where to send reports to, in addition to those in `config/emails.conf`.
Workflow describes the set of plugins, and their order.
Workflows are further described in [config/plugins.conf](PLUGINS.md).

## Optional

These files are found in the Illumina run directory although they are optional.

* CompletedJobInfo.xml
* config.xml
* QC/CompletedJobInfo.xml
* QC/GenerateFASTQRunStatistics.xml
* QC/RunInfo.xml
* QC/runParameters.xml
* QC/InterOp/ - this directory contains other files from the default Illumina run directory

