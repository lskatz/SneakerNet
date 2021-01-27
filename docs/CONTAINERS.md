# Containers

We have containerized SneakerNet into a docker image that is publicly available on [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).

Some test data can be found in this repository in
t/M00123-18-001-test/
and is generally what the variable `INDIR` would hold in this document.

## Requirements
Docker, Singularity,  or another Docker-compatible container software must be installed e.g. shifter (untested)

## Database(s)

As of SneakerNet version 0.11.0, the Kraken database must be installed separately.
You will need to choose a directory to keep the database in.
An example path is given in the instructions below, with `KRAKEN_DEFAULT_DB`.
```bash
KRAKEN_DEFAULT_DB=$HOME/db/kraken1/minikraken_20171013_4GB
mkdir -pv $KRAKEN_DEFAULT_DB
pushd $KRAKEN_DEFAULT_DB
cd ..
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz
tar -zxvf minikraken_20171019_4GB.tgz
rm -vf minikraken_20171019_4GB.tgz
popd
```

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
    singularity build sneakernet.sif docker://lskatz/sneakernet:latest

### Running SneakerNet using Singularity

Create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
The variable `MISEQ` represents the original MiSeq directory,
which is then converted into a new format in a new directory, `INDIR`.
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below, where file transfer and email is disabled.

    export MISEQ=my/miseq/run/dir
    export INDIR=my/run/dir

    # This line is here to catch you if you missed the KRAKEN_DEFAULT_DB step
    [[ -e "$KRAKEN_DEFAULT_DB" ]] || echo "ERROR: please set KRAKEN_DEFAULT_DB"
    
    # this assumes $MISEQ and $INDIR are in your $PWD
    singularity exec -B $PWD:/data sneakernet.sif SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
    singularity exec -B $PWD:/data -B $KRAKEN_DEFAULT_DB:/kraken-database sneakernet.sif SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR

## Docker

### Docker image installation

Docker CE must first be installed onto your system. Check that it is installed by running:
```bash
docker info
```
If it is not installed, visit: [https://docs.docker.com/install/](https://docs.docker.com/install/) and follow the install instructions according to your operating system.

For Ubuntu users, please see concise install instructions here: [https://github.com/StaPH-B/scripts/blob/master/image-information.md#docker-ce](https://github.com/StaPH-B/scripts/blob/master/image-information.md#docker-ce)

Download the docker image using the `latest` docker image tag
```bash
docker pull lskatz/sneakernet:latest
```

### Running SneakerNet using Docker

As in the Singularity section,
create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
The variable `MISEQ` represents the original MiSeq directory,
which is then converted into a new format in a new directory, `INDIR`.
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below, where file transfer and email is disabled.

```bash
export MISEQ=my/miseq/run/dir
export INDIR=my/sneakernet/run/dir

# This line is here to catch you if you missed the KRAKEN_DEFAULT_DB step
[[ -e "$KRAKEN_DEFAULT_DB" ]] || echo "ERROR: please set KRAKEN_DEFAULT_DB"

# -v flag will mount your PWD into the /data directory inside the container
# -u flag preserves your user/group when executing commands in the container
# --rm flag will remove/delete the container after it exits 
# make sure output files are written to /data so you don't lose them after the container exits!
docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
docker run --rm -v $PWD:/data -v $KRAKEN_DEFAULT_DB:/kraken-database -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR
```
