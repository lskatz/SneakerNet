# SneakerNet

A pipeline for processing reads from a sequencing run. Currently only for Illumina-based runs.

## Synopsis

What do you do with a MiSeq run after it finishes? Are there basic analyses that you run 
such as counting how many reads you obtained? Checking for contamination? **SneakerNet performs
all these initial analyses with a nice report that is emailed to you at the end.** Additionally,
there is a plugins directory such that each installation is expandable.

## Installation

See INSTALL.md

## Running it

A given run directory must have machine name, two-digit year, increment, and optionally a remark
in the name. These fields are delimited by a dash. For example, M1234-18-001-test.

Each run directory must also have a SampleSheet.csv file, fastq.gz files, and the InterOp directory.

### Testing for the impatient

If you have altered the configuration properly, then go ahead and test the software!

    $ SneakerNet.pl --test --now --numcpus 4

Also, please look at the `example` folder for another readme file and example dataset.

### Running it on an existing run directory in place

If you don't want to move the directory but just want to analyze it.  Let us say that
the test directory name is M1234-18-001-test.

    $ SneakerNetPlugins.pl --numcpus 4 M1234-18-001-test

