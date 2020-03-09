# SneakerNet Input folder

## Required

### samples.tsv

This file can be created from SampleSheet.csv by the plugin [sn_parseSampleSheet.pl](../docs/plugins/sn_parseSampleSheet.pl.md), from `SneakerNet.ro.pl --createsamplesheet`, or manually.
It has three columns, separated by tab: 

* sample - Name for a sample
* info - semicolon-delimited information. Keys are separated by an equals sign. For example: `taxon=Vibrio;route=calcEngine` means that the taxon is Vibrio and the ultimate destination for the files are in calcEngine. Vibrio is the key for the taxon in `config/taxonProperties.conf` and will not have meaning if there is no entry.
* filenames - semicolon-delimited fastq files. For example, `./2018AW-0585_S318_L001_R1_001.fastq.gz;./2018AW-0585_S318_L001_R2_001.fastq.gz`

## Optional

### snok.txt

### CompletedJobInfo.xml

### config.xml


### QC/CompletedJobInfo.xml

### QC/GenerateFASTQRunStatistics.xml

### QC/RunInfo.xml

### QC/runParameters.xml

### QC/InterOp/

