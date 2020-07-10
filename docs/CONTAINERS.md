# Containers

We have containerized SneakerNet into a docker image that is publicly available on [dockerhub](https://hub.docker.com/r/lskatz/sneakernet).

## TOC
- [Requirements](#requirements)
- [Singularity](#singularity)
  * [Singularity image installation](#singularity-image-installation)
  * [Running SneakerNet using Singularity](#running-sneakernet-using-singularity)
- [Docker](#docker)
  * [Docker image installation](#docker-image-installation)
  * [Running SneakerNet using Docker](#running-sneakernet-using-docker)

## Requirements
Docker, Singularity,  or another Docker-compatible container software must be installed e.g. shifter (untested)

## Singularity

### Singularity image installation
Check to make sure Singularity is installed:
```bash
singularity --help
```
If it is not installed, please visit the following link for installing the latest version of Singularity (at time of writing this). You will need administrator priveleges and to install some additional dependencies before installing singularity itself. [https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps](https://sylabs.io/guides/3.4/user-guide/quick_start.html#quick-installation-steps)


(OPTIONAL) Make a fast scratch directory if your current drive is not fast.
This step, including the `SINGULARITY_TMPDIR` variable is optional.

    export SINGULARITY_TMPDIR=/scratch/$USER/singularity-tmp
    mkdir -pv $SINGULARITY_TMPDIR

Navigate into a directory where you would like your Singularity image to be stored.
Build the image with `singularity build`.

    cd your/containers/directory/
    singularity build sneakernet.simg docker://lskatz/sneakernet:latest

### Running SneakerNet using Singularity

Create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below where file transfer and email is disabled.
Commands are broken into multiple lines for readability.
```bash
# Set up a SneakerNet style directory using the example data
singularity exec -B $(pwd):/data sneakernet.simg \
SneakerNet.roRun.pl /SneakerNet-*/t/M00123-18-001-test -o /data/singularity-test

# Run SneakerNet on the example data
singularity exec -B $(pwd):/data sneakernet.simg \
SneakerNetPlugins.pl --numcpus 8 --no email --no transfer --no save /data/singularity-test

#####################################

# Run SneakerNet on your own data (typically a MiSeq run directory)

# this assumes INDIR and OUTDIR are in your $PWD
export INDIR=my-input-miseq-run-dir/
export OUTDIR=my-output-dir/

# Set up a SneakerNet style directory using your own data
singularity exec -B $(pwd):/data sneakernet.simg \
SneakerNet.roRun.pl /data/$INDIR -o /data/$OUTDIR

# Run SneakerNet on your own data
singularity exec -B $(pwd):/data sneakernet.simg \
SneakerNetPlugins.pl --numcpus 8 --no email --no transfer --no save /data/$OUTDIR
```

## Docker

### Docker image installation

Docker CE must first be installed onto your computer. Check that it is installed by running:
```bash
docker info
```
If it is not installed, visit: [https://docs.docker.com/install/](https://docs.docker.com/install/) and follow the install instructions according to your operating system.

Download the docker image using the `latest` docker image tag
```bash
docker pull lskatz/sneakernet:latest
```

### Running SneakerNet using Docker

Commands broken into multiple lines for readability.
```bash
# Docker run options explanation:
# docker run -v option will mount your PWD into the /data directory inside the container
# docker run -u option preserves your user/group when executing commands in the container
# docker run --rm option will remove/delete the container after it exits 

# Set up a SneakerNet style directory using the example data
# Make sure output files are written to /data so you don't lose them after the container exits!
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data lskatz/sneakernet:latest \
SneakerNet.roRun.pl /SneakerNet-*/t/M00123-18-001-test -o /data/test

# Run SneakerNet on example data
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data lskatz/sneakernet:latest \
SneakerNetPlugins.pl --numcpus 8 --no email --no transfer --no save /data/test

#####################################

# Run SneakerNet on your own data (typically a MiSeq run directory)

# this assumes INDIR and OUTDIR are in your $PWD
export INDIR=my-input-miseq-run-dir/
export OUTDIR=my-output-dir/

# Set up a SneakerNet style directory using your own data
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data lskatz/sneakernet:latest \
SneakerNet.roRun.pl /data/$INDIR -o /data/$OUTDIR

# Run SneakerNet on your own data
docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data lskatz/sneakernet:latest \
SneakerNetPlugins.pl --numcpus 8 --no email --no transfer --no save /data/$OUTDIR
```
