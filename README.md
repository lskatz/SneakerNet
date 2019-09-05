# SneakerNet

A pipeline for processing reads from a sequencing run. Currently only for Illumina-based runs.

## Synopsis

What do you do with a MiSeq run after it finishes? Are there basic analyses that you run 
such as counting how many reads you obtained? Checking for contamination? **SneakerNet performs
all these initial analyses with a nice report that is emailed to you at the end.** Additionally,
there is a plugins directory such that each installation is expandable.

## Installation

See INSTALL.md

## Workflow

### Creating a SneakerNet project directory

SneakerNet requires a project directory that is in a certain format already.
To create the project, you can use `SneakerNet.ro.pl`.  For example,

    SneakerNet.ro.pl --createsamplesheet -o M1234-18-001-test miseq/working/directory

M01234-19-01-test is a project folder name, where it is dash-delimited and contains
machine name, year, ordinal, and optionally a name.

### Running SneakerNet

Because SneakerNet takes a long time to run through all plugins, it is
a good idea to pipe the output to a file and then follow it with `tail -f`.

    SneakerNetPlugins.pl --numcpus 8 M1234-18-001-test > M1234-18-001-test/SneakerNet.log 2>&1 &
    tail -f M1234-18-001-test/SneakerNet.log

## Output

SneakerNet produces a subfolder `SneakerNet` in your run directory.
It also emails a report. To view a sample report, please go to example/M00123-18-001-test/SneakerNet/forEmail/report.html 
in this repository.

## Plugins

SneakerNet is based on plugins.  In this context, a plugin is an independent script
that can run an analysis on a run directory, accept standard inputs (e.g., `--help`),
and create standard output files.

### Plugins for developers

Please look at the (readme for plugins)[SneakerNet.plugins/README.md].

