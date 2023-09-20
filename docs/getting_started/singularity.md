# Running Labelmerge with Singularity

## Pre-requisites
1. Singularity / Apptainer is is installed on your system. For more info, see
the detailed [Apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
1. The following command-line tools are installed:
    * wget
1. Sufficient disk-space (rough estimate)
1. Sufficient CPU and memory - We recommend at least 16GB memory if using default parameters.

## First time setup
Pull the container. This can be done from DockerHub, but requires a large 
amount of disk space in your `/tmp` folder, since it has to convert from a 
Docker container to a Singularity/Apptainer container. The example below pulls
the latest versioned container (replace `latest` with `vX.X.X` for a specific
version).

```
singularity pull docker://khanlab/labelmerge:latest
```
_Note: If you encounter any errors pulling the container from DockerHub, it may
be because you are running out of disk space in your cache folders. You can 
change these locations by setting environment variables, however, using a 
network file system for the folders may result in poor performance:_
```
export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity
```


Run Labelmerge without any arguments to print the short help:

```
singularity run -e khanlab_labelmerge_latest.sif
```

Use the `-h` option to get a detailed help listing:

```
singularity run -e khanlab_labelmerge_latest.sif -h
```

_Note that all the Snakemake command-line options are also available,
and can be listed with `--help-snakemake`:

```
singularity run -e khanlab_labelmerge_latest.sif --help-snakemake
```

### Explanation

Everything prior to the container (`khanlab_labelmerge_latest.sif`) are arguments
to Singluarity / Apptainer, and after are to Labelmerge itself. The first three arguments to 
Labelmerge (as with any BIDS App) are the input folder (`test/data/bids`), the 
output folder (`test/data/derivatives`), and the analysis level (`participant`).
The `participant` analysis level is used in Labelmerge to perform further 
participant-level processing of external atlases (combining segmentations). 

```
singularity run -e khanlab_labelmerge_latest.sif test/data/bids test/data/derivatives participant -np
```

Now to actually run the workflow, we need to specify how many cores to use and 
leave out the dry-run option. The Snakemake `--cores` option tells Labelmerge how
many cores to use. Using `--cores 8` means that Labelmerge will only make use of 8 
cores at most. Generally speaking, you should use `--cores all`, so it can make 
maximal use of all available CPU cores it has access to on your system. This is 
especially useful if you are running multiple subjects.


_Note that you may need to adjust your 
[Singularity / Apptainer options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) 
to ensure the container can read and write to your input and output directories, 
respectively. You can bind paths easily by setting an environment variable, 
e.g. if you have a `/project` folder that contains your data, you can add it to
the `SINGULARITY_BINDPATH` so it is available when you are running a container:_

```
export SINGULARITY_BINDPATH=/data:/data
```

After this completes, you have additional folders in your output folder,
`test/data/derivatives`, for the one subject.
