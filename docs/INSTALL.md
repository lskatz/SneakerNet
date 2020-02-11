# Installation

## Quick installation

    mkdir ~/bin
    cd bin
    git clone https://github.com/lskatz/SneakerNet.git
    cd SneakerNet
    make

## Test the installation

    make test

## Dependencies

* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* Multithreaded Perl (already installed on most computers)
* Kraken: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* _To be continued..._

## Configuration

You will need to edit some files for configuration before using SneakerNet.

    $ cp -r config.bak config
    $ cd config

### emails

List any emails here, comma-separated. These emails will be sent reports by default for each
SneakerNet run.

### taxonProperties

Each taxon that you might have sequences for is defined here. If not defined here, nothing bad
will happen though.  For each taxon, you can show the minimum quality and coverage thresholds;
the genome size in bp; the regular expression to match filename to taxon; and the destination
folder on the remote computer.

### settings

This file has certain key/values and should be left alone, unless you are a developer.

### plugins

This file defines workflows.
The first keyword is the workflow. After the equals sign,
a number of plugins are listed and delimited by commas.
They are shown in order of execution. For example, all
assembly-based plugins must appear after the genome
assembly plugin runs.

You are able to create your own workflow or use the 
pre-defined ones here.

Each workflow is listed in the SneakerNet `snok.txt` file
in the format of `workflow = default` where in this example
the default workflow is invoked. If no workflow is given
in `snok.txt` or if `snok.txt` is missing, the default
workflow will be used.

