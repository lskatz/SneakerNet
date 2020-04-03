# Containers

We have containerized SneakerNet into [dockerhub](https://hub.docker.com/repository/docker/lskatz/sneakernet).

## Singularity

### Singularity installation

Make a fast scratch directory if your current drive is not fast.
This step, including the `SINGULARITY_TMPDIR` variable is optional.

    export SINGULARITY_TMPDIR=/scratch/$USER/singularity-tmp
    mkdir -pv SINGULARITY_TMPDIR

Navigate into a directory where you would like your Singularity file.
Build the image with `singularity build`.

    cd your/containers/directory/
    singularity build sneakernet.simg docker://lskatz/sneakernet:latest

### Singularity running

Create your SneakerNet-formatted directory first.
For instructions on how to create the SneakerNet-formatted directory, please see [SneakerNetInput.md](SneakerNetInput.md).
Next, run `SneakerNetPlugins.pl` on the target directory.
An example is below where file transfer and email is disabled.

    export INDIR=my/run/dir
    SneakerNet.roRun.pl my/MiSeq/run/directory -o $INDIR
    singularity exec -B $INDIR:/data ../../containers/sneakernet.simg SneakerNetPlugins.pl --numcpus 12 --no email --no transfer --no save $INDIR

## Docker

