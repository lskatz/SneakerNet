# SneakerNet

## Synopsis

A pipeline for processing reads from a sequencing run. Currently supports Illumina or Ion Torrent,
but it can be expanded to other platforms.

    # Run SneakerNet on the example data
    SneakerNetPlugins.pl --numcpus 4 t/data/M00123-18-001

<p align='center'>
  <img src='./docs/images/overview.png' alt='SneakerNet workflow' width='400' />
</p>

## Installation

See [docs/INSTALL.md](docs/INSTALL.md)

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
It also emails a report. To view a sample report, please go to 
[example/M00123-18-001-test/SneakerNet/forEmail/report.html]
in this repository.

## Plugins

SneakerNet is based on plugins.  In this context, a plugin is an independent script
that can run an analysis on a run directory, accept standard inputs (e.g., `--help`),
and create standard output files.

For more details, see the [plugins readme](docs/PLUGINS.md).

### Plugins for developers

Please look at the [readme for plugins](docs/PLUGINSDEV.md).

