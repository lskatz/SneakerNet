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

There are few core dependencies. However, SneakerNet is
based on plugins that have individual dependencies.
Here is a list of dependencies across the board, although
it is possible that all of the following are not needed
for a complete SneakerNet analysis.

To find the right dependencies for you, run the following

    SneakerNet.checkdeps.pl --list # and find the workflow, e.g., 'default'
    SneakerNet.checkdeps.pl default

### Comprehensive list

This list was created using `SneakerNet.checkdeps.pl iontorrent metagenomics cryptosporidium default`

* Multithreaded Perl (already installed on most computers)
* GNU utilities (`cp`, `cat`, ...)
* `rsync`
* `ssh`
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* Skesa
* Prodigal
* Shovill
* `mlst`
* ColorID
* `zip`

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

