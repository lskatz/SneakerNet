# Containers

We have containerized SneakerNet into a docker image that is publicly available on [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).

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

    export MISEQ=my/miseq/run/dir
    export INDIR=my/run/dir
    
    # this assumes $MISEQ and $INDIR are in your $PWD
    singularity exec -B $PWD:/data sneakernet.simg SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
    singularity exec -B $PWD:/data sneakernet.simg SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR

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

Run the container, using the same example as above:
```bash
# this assumes that your miseq run directory, and desired SneakerNet input directory are in your PWD.
# you will need to adjust PATHs depending on where your data is located
export MISEQ=my/miseq/run/dir
export INDIR=my/sneakernet/run/dir

# -v flag will mount your PWD into the /data directory inside the container
# -u flag preserves your user/group when executing commands in the container
# --rm flag will remove/delete the container after it exits 
# make sure output files are written to /data so you don't lose them after the container exits!
docker run --rm -v $PWD:/data -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNet.roRun.pl /data/$MISEQ -o /data/$INDIR
docker run --rm -v $PWD:/data -u $(id -u):$(id -g) lskatz/sneakernet:latest SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save /data/$INDIR
```
