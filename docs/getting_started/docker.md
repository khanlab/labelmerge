# Running Labelmerge with Docker on Windows

**Note, these instructions assume you have Docker already installed on a Windows system.
Docker can also be run on Linux or MacOS with similar commands, but here, we 
will assume the default Windows CLI is being used.

## First time setup

Open your Windows Command Prompt by clicking the `Windows` button and typing
`cmd` and pressing the `Enter` on your keyboard. This is where you will enter 
your SCATTR commands. Feel free to make a new directory with `mkdir` or move to
a directory you would like to work out of with `cd. For this example, we will
work from:

```
cd c:\Users\username\Downloads
```

Pull the container (this will take some time and storage stage, but like an 
installation, it only needs to be done once and can be then be run on many 
datasets). The example below pulls the latest versioned container (replace 
`latest` with `vX.X.X` for a specific version).

```
docker pull khanlab/labelmerge:latest
```

Run without any arguments to print hte short help:

```
docker run -it --rm khanlab/labelmerge:latest
```

Use the `-h` option to get a detailed help listing:

```
docker run -it --rm khanlab_labelmerge_latest.sif -h
```

_Note that all the Snakemake command-line options are also available in SCATTR,
and can be listed with `--help-snakemake`:

```
docker run -it --rm khanlab_labelmerge_latest.sif --help-snakemake
```

## Running an example

** - Dummy Example for Scattr - We will use the `test` folder found from the 
[Github repository](https://github.com/khanlab/scattr/tree/main/test/) via
`git clone` to the previously mentioned folder to demonstrate an example of 
how to run SCATTR

```
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test khanlab_scattr_latest.sif /test/data/bids /test/data/derivatives participant --fs-license /test/fs_license --force-output -n
```

### Explanation

Everything prior to the container (`khanlab_labelmerge_latest.sif`) are arguments
to Docker and after are to Lablemerge itself. The first three arguments to Docker
are to enable interactive mode (`-it`), run and subsequently remove the Docker
container upon completion (`--rm`) and mount the the directoty 
(`-v c:\Users\username\Downloads\labelmerge\test`) to a directory within the
container named `\test`. These are not specific to Labelmerge, but are general ways
to use Docker. You may want to familiarize yourself with 
[Docker options](https://docs.docker.com/engine/reference/run/).

The first three arguments to Labelmerge (as with any BIDS App) are the input folder 
(`test/data/bids`), the output folder (`test/data/derivatives`), and the 
analysis level (`participant`). The `participant` analysis level is used in 
SCATTR to perform further participant-level processing of Freesurfer (thalamus 
segmentation), external atlases (combining segmentations) and diffusion derived
data (estimation of fibre orientation distributions). This includes estimating 
an average response function from input data. The `--fs-license` argument allows
for specification of the location of the required Freesurfer license file if not 
already specified in the `FS_LICENSE` environment variable. Note, that 
this is required to perform any Freesurfer-related processing. The 
`--force-output` flag is a Snakemake argument that is invoked to allow for 
writing of output file to already existing folders - in this case, for thalamus 
segmentations via Freesurfer. We also used the `--dry-run/-n` option to print 
out what would run, without actually running the workflow.

When you run the above command, a long listing will print out, describing all 
the rules that will be run. We can also have a shell command used for each rule
printed to screen using the `-p` Snakemake option

```
docker run -it --rm -v c:\Users\username\Downloads\scattr\test:\test  khanlab_scattr_latest.sif /test/data/bids /test/data/derivatives participant --fs-license /test/fs_license --force-output -np
```

Now to actually run the workflow, we need to specify how many cores to use and 
leave out the dry-run option. The Snakemake `--cores` option tells SCATTR how
many cores to use. Using `--cores 8` means that SCATTR will only make use of 8 
cores at most. Generally speaking, you should use `--cores all`, so it can make 
maximal use of all available CPU cores it has access to on your system. This is 
especially useful if you are running multiple subjects.

After this completes, you have additional folders in your output folder,
`c:\Users\username\Downloads\labelmerge\test\data\derivatives`, for the one subject.

### Exploring different options

