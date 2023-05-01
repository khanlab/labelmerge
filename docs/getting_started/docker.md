# Running Labelmerge with Docker on Windows

**Note, these instructions assume you have Docker already installed on a Windows system.
Docker can also be run on Linux or MacOS with similar commands, but here, we 
will assume the default Windows CLI is being used.

## First time setup

Open your Windows Command Prompt by clicking the `Windows` button and typing
`cmd` and pressing the `Enter` on your keyboard. This is where you will enter 
your SCATTR commands. Feel free to make a new directory with `mkdir` or move to
a directory you would like to work out of with `cd'. For this example, we will
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

Run without any arguments to print the short help:

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
