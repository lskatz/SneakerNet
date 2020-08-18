# Installation

## Quick installation

    mkdir ~/bin
    cd bin
    git clone https://github.com/lskatz/SneakerNet.git
    cd SneakerNet
    
    # the following two lines are for local installations
    cpanm --local-lib=$HOME local::lib && eval $(perl -I $HOME/lib/perl5/ -Mlocal::lib)
    export PERL5LIB=$PERL5LIB:$HOME/lib/perl5:$HOME/lib/perl5/x86_64-linux-gnu-thread-multi:$HOME/lib/perl5/x86_64-linux-gnu-thread-multi/auto
    
    # The following lines are regardless of local or global installation
    cpanm version::vpp
    cpanm --installdeps --notest . -l .
    perl Makefile.PL
    make

## Test the installation

    make test

## Dependencies

## Resources
The minimum requirements are about 4G RAM, but 8G RAM is recommended.
An average run with only 1 thread is about 200 minutes; however with 4 threads, the run time gets down to about 9 minutes.
For more information: https://github.com/lskatz/SneakerNet/issues/42#issuecomment-672900992

## Software
There are few core dependencies. However, SneakerNet is
based on plugins that have individual dependencies.
Here is a list of dependencies across the board, although
it is possible that all of the following are not needed
for a complete SneakerNet analysis.

To find the right dependencies for you, run the following

    SneakerNet.checkdeps.pl --list  # and find the workflow, e.g., 'default'
    SneakerNet.checkdeps.pl default # Check deps for all plugins for the workflow 'default'

#### Comprehensive list

This list was created using `SneakerNet.checkdeps.pl [iontorrent metagenomics cryptosporidium default]`
on version 0.8.14.
Dependencies may or may not have changed since then but they can still be checked using this script.

All workflows require perl v5.12 or higher, compiled with multithreading.
This version of perl is already installed in most modern Linux operating systems.

##### Default workflow

* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* Python3
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* sendmail
* Skesa
* staramr
* `zip`

##### metagenomics

* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* `zip`

##### cryptosporidium

* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* Shovill
* `countGP60repeats.pl`: currently in development in a private repo. To exclude, remove `sn_crypto_gp60.pl` from `config/plugins.conf` (already not included by default)
* `zip`

##### iontorrent

* Shovill
* SPAdes
* blastn (BLAST+)
* GNU utilities (`cp`, `cat`, ...)
* Kraken1: http://ccb.jhu.edu/software/kraken/
* Krona: https://github.com/marbl/Krona/
* `mlst`: https://github.com/tseemann/mlst
* Prodigal
* Python3
* CG-Pipeline: https://github.com/lskatz/cg-pipeline (scripts only -- do not run `make`)
* sendmail
* Skesa
* staramr
* `zip`

##### Non-workflow

###### SneakerNet.roRun.pl

* bcl2fastq

## Configuration

You will need to edit some files for configuration before using SneakerNet.

    $ cp -r config.bak config
    $ cd config

Some settings are necessary to change in the config/\*.conf files.
Please edit these files accordingly.

### emails

List any emails here, comma-separated. These emails will be sent reports by default for each
SneakerNet run.

### taxonProperties

_see [INSTALL.taxonProperties.md](INSTALL.taxonProperties.md) for more information_

Each taxon that you might have sequences for is defined here. If not defined here, nothing bad
will happen though.  For each taxon, you can show the minimum quality and coverage thresholds;
the genome size in bp; the regular expression to match filename to taxon; and the destination
folder on the remote computer.

### settings

This file has certain key/values and should be carefully edited.
Some individual plugin settings are found here.

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

