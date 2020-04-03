# SYNOPSIS

`sn_immediateStatus.pl`

Emails an immediate report on whether a run has warnings

# Software requirements

# Algorithm

This plugin is meant to be run as one of the first plugins in a workflow.
It reads `samples.tsv` and the fastq files in the directory and determines:

1. Does every sample have fastq file(s)?
2. Is every fastq file associated with a sample?

# Outputs

## Table

Table of the immediate reaction, with columns

* ErrType
  * fastq - A fastq file was found, but no sample "owns" it.
  * sample - A sample was listed in `samples.tsv` but no fastq files were found.
* Sample - either the sample or the fastq file in question
* ErrKeyword - A shorthand way of writing out the Error column
* Error - A specific error message. It might offer help such as suggesting the correct fastq files or the correct sample.

## email

An email with the table is sent to those listed in `snok.txt` and `emails.conf`.

