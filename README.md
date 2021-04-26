# SneakerNet

## Synopsis

A pipeline for processing reads from a sequencing run. Currently supports Illumina or Ion Torrent,
but it can be expanded to other platforms.

    # Run SneakerNet on the example data
    SneakerNetPlugins.pl --numcpus 4 t/M00123-18-001-test

<p align='center'>
  <img src='./docs/images/overview.png' alt='SneakerNet workflow' width='400' />
</p>

### Main steps

This is the default workflow in v0.14
but there are other workflows available as described in
[PLUGINS.md](/docs/PLUGINS.md#workflows).
 
* [Parse sample entries](/docs/plugins/sn_detectContamination-kraken.pl.md) - create an input file `samples.tsv`
* [Read metrics](/docs/plugins/addReadMetrics.pl.md) - get raw read yields and raw read QC summary (CG-Pipeline)
* [Assembly](/docs/plugins/assembleAll.pl.md) - assemble each genome (Shovill/skesa)
* [MLST](/docs/plugins/sn_mlst.pl.md) - 7-gene MLST (_mlst_)
* [Run Kraken](/docs/plugins/sn_kraken.pl.md)
* [Contamination detection](/docs/plugins/sn_detectContamination-kraken.pl.md) - check that all reads come from one taxon for each genome (Kraken)
* [Contamination detection](/docs/plugins/sn_detectContamination-mlst.pl.md) - check that all seven MLST genes have only one instance in the genome as expected (ColorID)
* [Base balance](/docs/plugins/baseBalance.pl.md) - check that the ratio of A/T is approximately 1 and same with C/T
* [Antimicrobial resistance gene prediction](/docs/plugins/sn_staramr.pl.md) - detect genotype and predict phenotype (staramr)
* [Pass/fail](/docs/plugins/sn_passfail.pl.md) - list all genomes that have failed Q/C
* [Transfer Files](/docs/plugins/transferFilesToRemoteComputers.pl.md) - files are copied to a remote folder
* [HTML summary report](/docs/plugins/sn_report.pl.md)
* [Email](/docs/plugins/emailWhoever.pl.md) the report

### Quick start

1. Install and configure SneakerNet - [from source](docs/INSTALL.md) or [with a container](docs/CONTAINERS.md)
2. Make an input folder from your MiSeq run [docs/SneakerNetInput.md](docs/SneakerNetInput.md)
3. Run [`SneakerNetPlugins.pl`](docs/SneakerNetPlugins.pl.md) on the input folder.

## Installation

See [docs/INSTALL.md](docs/INSTALL.md)

_NOTE_: to ensure all dependencies are met, please follow
the dependencies section under the [installation document](docs/INSTALL.md).

### Container installation

SneakerNet has been containerized and is at [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).
For more information, please see our [containers documentation](docs/CONTAINERS.md).

Here is a summary of Docker commands, from the [containers documentation](docs/CONTAINERS.md).

    # Pull image
    docker pull lskatz/sneakernet:latest
    # Import data directly from the MiSeq machine, where $MISEQ is a raw run folder exported by the MiSeq machine
    # and $INDIR is the newly created SneakerNet input folder
    docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
    # Run SneakerNet on the $INDIR (SneakerNet formatted folder)
    docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR

## Workflow

### Creating a SneakerNet project directory

_For more information on a SneakerNet-style folder, see [docs/SneakerNetInput.md](docs/SneakerNetInput.md)_

SneakerNet requires a project directory that is in a certain format already.
To create the project, you can use `SneakerNet.roRun.pl`.  For example,

    SneakerNet.roRun.pl --createsamplesheet -o M1234-18-001-test miseq/working/directory

M01234-19-01-test is a project folder name, where it is dash-delimited and contains
machine name, year, ordinal, and optionally a name.
Fastq files must be in the format of `_R1_` instead of `_1` and `_R2_` instead of `_2` for this particular script to parse the files properly.

### Running SneakerNet

It is generally a good idea to edit a file `snok.txt` to configure the run further.
For more information on the workflow, see the configuration section in `INSTALL.md`.
For example,

    echo "emails = example@example.com, blah@example.com" > t/data/M00123-18-001/snok.txt
    echo "workflow = default" >> t/data/M00123-18-001/snok.txt

And then run SneakerNet like so (optionally following the log with `tail -f`):

    SneakerNetPlugins.pl --numcpus 8 t/data/M00123-18-001 > t/data/M00123-18-001/SneakerNet.log 2>&1 &
    tail -f t/data/M00123-18-001/SneakerNet.log

#### Containers

SneakerNet has been containerized and is at [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).
For more information, please see our [containers documentation](docs/CONTAINERS.md).

## Output

_For more information, please see [docs/SneakerNetOutput.md](docs/SneakerNetOutput.md)_

SneakerNet produces a subfolder `SneakerNet/` in your run directory.
It also emails a report. To view a sample report, please go to 
t/report.html
in this repository.

## Plugins

SneakerNet is based on plugins.  In this context, a plugin is an independent script
that can run an analysis on a run directory, accept standard inputs (e.g., `--help`),
and create standard output files.

For more details, see the [plugins readme](docs/PLUGINS.md).

### Plugins for developers

You too can develop for SneakerNet!  For more information, 
please look at the [readme for plugins](docs/PLUGINSDEV.md)
and the [contributing](CONTRIBUTING.md) doc.

## Further reading

Please see the docs subfolder for more specific documentation.

