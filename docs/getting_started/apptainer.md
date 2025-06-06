# Running Labelmerge with Singularity

## Pre-requisites
1. Apptainer (or Singularity) is installed on your system. For more info, see the detailed [apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
 2. The following command-line tools are installed:
      - wget
      - tar
 3. Sufficient disk-space is needed
 4. Sufficient CPU and memory - the more you have, the faster it will run, but we recommend at least 4 CPU cores and 16GB memory.


## First time setup
Pull the container. This can be done from DockerHub, but requires a large 
amount of disk space in your `/tmp` folder, since it has to convert from a 
Docker container to a Singularity/Apptainer container. The example below pulls
the latest versioned container (replace `latest` with `vX.X.X` for a specific
version).

Pull the container:
```
apptainer pull khanlab_labelmerge_latest.sif docker://khanlab/labelmerge:latest
```

Run Labelmerge without any arguments to print the short help:

```
apptainer run -e khanlab_labelmerge_latest.sif
```

Use the `-h` option to get a detailed help listing:

```
apptainer run -e khanlab_labelmerge_latest.sif -h
```

_Note that all the Snakemake command-line options are also available,
and can be listed with `--help-snakemake`:

```
apptainer run -e khanlab_labelmerge_latest.sif --help-snakemake
```

Note: If you encounter any errors pulling the container from dockerhub, it may be because you are running 
out of disk space in your cache folders. Note, you can change these locations 
by setting environment variables, however, using a network file system for the folders may result in poor performance and/or errors e.g.:
    
```
export APPTAINER_CACHEDIR=/YOURDIR/.cache/apptainer
```

## Running an example

You can try Labelmerge on a sample dataset to make sure everything works as expected.

First, download and extract a single-subject BIDS dataset for this test:

```bash
wget "https://www.dropbox.com/scl/fo/qzsym6f7k56yc8jcseawu/AOiV4AhmH6oiTO0Cr4GoaW8?rlkey=wftu1ph2cbdlysocqvbn1muka&st=xw5x84bp&dl=0" -O labelmerge_test.zip
unzip labelmerge_test.zip
```

This comand will create a `/tpl-MNI152NLin2009cAsym` folder with data from the **MNI152NLin2009cAsym_atlas**, containing cortical and subcortical volumes to merge.

```
tpl-MNI152NLin2009cAsym
    └── anat
           ├── tpl-MNI152NLin2009cAsym_atlas-MIAL67ThalamicNuclei_desc-tn_dseg.tsv
           ├── tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.tsv
           ├── tpl-MNI152NLin2009cAsym_res-01_atlas-MIAL67ThalamicNuclei_desc-tn_dseg.nii.gz
           └──  tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.nii.gz

3 directories, 4 files
```

Now let's run labelmerge. 

    apptainer run -e khanlab_labelmerge_latest.sif tpl-MNI152NLin2009cAsym/ output_dir participant --base-desc 100Parcels7Networks --overlay_bids_dir  tpl-MNI152NLin2009cAsym/ --overlay_desc tn -n

### Explanation
Everything prior to the container (`khanlab_labelmerge_latest.sif`) are arguments to apptainer, and after are to labelmerge itself. The first three arguments to labelmerge (as with any BIDS App) are the input
folder (`/tpl-MNI152NLin2009cAsym`), the output folder (`output_dir`), and then the analysis level (`participant`). The `participant` analysis 
level is used in labelmerge for performing any
participant-level processing. We also used the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.

When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    apptainer run -e khanlab_labelmerge_latest.sif  labelmerge_test output_dir participant --base-desc 100Parcels7Networks --overlay_bids_dir  labelmerge_test --overlay_desc tn -np | less

Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells labelmerge how many cores to use.
 Using `--cores 8` means that labelmerge will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

    apptainer run -e khanlab_labelmerge_latest.sif  labelmerge_test output_dir participant --base-desc 100Parcels7Networks --overlay_bids_dir  labelmerge_test --overlay_desc tn --cores all

Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/apptainer_run.html) to ensure the container can read and write to yout input and output directories, respectively. You can bind paths easily by setting an 
environment variable, e.g. if you have a `/project` folder that contains your data, you can add it to the `APPTAINER_BINDPATH` so it is available when you are running a container:

```
    export APPTAINER_BINDPATH=/data:/data
```

After this completes, you should have a `output_dir` folder with outputs for the one subject.
