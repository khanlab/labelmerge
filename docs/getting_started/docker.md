# Running Labelmerge with Docker on Windows

**Note, these instructions assume you have Docker already installed on a Windows system.
Docker can also be run on Linux or MacOS with similar commands, but here, we 
will assume the default Windows CLI is being used.**

## First time setup

Open your Windows Command Prompt by clicking the `Windows` button and typing
`cmd` and pressing the `Enter` on your keyboard. This is where you will enter 
your commands. Feel free to make a new directory with `mkdir` or move to
a directory you would like to work out of with `cd'. For this example, we will
work from:

```
cd c:\Users\username\Downloads
```

Pull the container (this will take some time and storage stage, but like an 
installation, it only needs to be done once and can then be run on many 
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

*Note that all the Snakemake command-line options are also available,
and can be listed with `--help-snakemake`:*

```
docker run -it --rm khanlab_labelmerge_latest.sif --help-snakemake
```

## Running an example

Download and extract a BIDS dataset for this test from [labelmerge_test.tar]("https://www.dropbox.com/scl/fo/qzsym6f7k56yc8jcseawu/AOiV4AhmH6oiTO0Cr4GoaW8?rlkey=wftu1ph2cbdlysocqvbn1muka&st=xw5x84bp&dl=0"). Here we will also assume you chose to save and extract to the directory `c:\Users\msnyder\Downloads\`.

This contains a `/labelmerge_test` directory with data from the **MNI152NLin2009cAsym_atlas**, containing cortical and subcortical volumes to merge.

```
labelmerge-test
    └── tpl-MNI152NLin2009cAsym
        └── anat
            ├── tpl-MNI152NLin2009cAsym_atlas-MIAL67ThalamicNuclei_desc-tn_dseg.tsv
            ├── tpl-MNI152NLin2009cAsym_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.tsv
            ├── tpl-MNI152NLin2009cAsym_res-01_atlas-MIAL67ThalamicNuclei_desc-tn_dseg.nii.gz
            └──  tpl-MNI152NLin2009cAsym_res-01_atlas-Schaefer2018_desc-100Parcels7Networks_dseg.nii.gz

3 directories, 4 files
```

Now let's run labelmerge on the test dataset. Docker will need read/write access to the input and output directories, respectively. This is achieved with the `-v` flag. This 'binds' or 'mounts' a directory to a new directory inside the container.

    docker run -it --rm -v c:\Users\msnyder\Downloads\labelmerge_test:/bids -v c:\Users\msnyder\Downloads\labelmerge_test_output:/output khanlab/labelmerge:latest /bids /output participant --base-desc 100Parcels7Networks --overlay_bids_dir /bids --overlay_desc tn -n

### Explanation

-v c:\Users\msnyder\Downloads\labelmerge_test:/bids` tells Docker to mount the directory `c:\Users\msnyder\Downloads\labelmerge_test` into a new directory inside the container named `/bids`. We then do the same things for our output directory named `labelmerge_test_output`, which we mount to `/output` inside the container. These arguments are not specific to labelmerge but rather are general ways to use Docker. You may want to familiarize yourself with [Docker options](https://docs.docker.com/engine/reference/run/).

Everything after we specified the container (`khanlab/labelmerge:latest`) are arguments to AutoAFIDs itself. The first of these arguments (as with any BIDS App) are the input directory (`/bids`), the output directory (`/output`), and then the analysis level (`participant`). The `participant` analysis 
level is used in labelmerge for performing any participant-level processing. We then need to specify the description on the base image and tsv for labelmerge with `--base-desc`, the overlay nifti bids directory with `--overlay_bids_dir` and the overlay image and tsv description with `--overlay-desc`. We also used the `--dry-run/-n`  option to just print out what would run, without actually running anything.

When you run the above command, a long listing will print out, describing all the rules that 
will be run. Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells labelmerge how many cores to use.
 Using `--cores 8` means that labelmerge will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

    docker run -it --rm -v c:\Users\msnyder\Downloads\labelmerge_test:/bids -v c:\Users\msnyder\Downloads\labelmerge_test_output:/output khanlab/labelmerge:latest /bids /output participant --base-desc 100Parcels7Networks --overlay_bids_dir tpl-MNI152NLin2009cAsym --overlay_desc tn --cores all

After this completes, you have a labelmerge_test_output directory with outputs for one subject. 

