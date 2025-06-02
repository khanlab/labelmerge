# Running Labelmerge with Conda

Labelmerge can be installed and run using Conda on **Linux and macOS** systems.

**Note:** Conda installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.

---

## For Users: Installing Labelmerge via Conda

These steps are intended for **end users** who simply want to run Labelmerge.

### 1. Install Conda (if not already installed)

Follow the instructions at the official Conda documentation:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

---

### 2. Create and activate a new Conda environment

```bash
conda install mamba -c conda-forge
mamba create --name labelmerge-env -c khanlab -c conda-forge -c bioconda labelmerge
eval "$(mamba shell hook --shell zsh)"
mamba activate labelmerge-env
```

---

### 3. Test the installation

Run the following command to verify the installation:

```bash
labelmerge -h
```

You should see a help message listing all available command-line options.  
If this runs successfully, you’re ready to start processing data with Labelmerge!

## Running an example

You can try Labelmerge on a sample dataset to make sure everything works as expected.

First, download and extract a single-subject BIDS dataset for this test:

***need to add test data to dropbox***

### Run the full Labelmerge BIDS pipeline

By default (Linux or Intel-based macOS), you can run:

```bash
labelmerge bids_dir output_dir participant --base-desc base_desc --overlay_bids_dir overlay_bids_dir --overlay_desc overlay_desc --cores all
```

This should run the full pipeline and place results in a new `output_dir/` folder.

If you’re on an M-chip mac, prefix with CONDA_SUBDIR=osx-64 to ensure compatibility:

```bash
CONDA_SUBDIR=osx-64 labelmerge bids_dir output_dir participant --base-desc base_desc --overlay_bids_dir overlay_bids_dir --overlay_desc overlay_desc --cores all
```

## Cache Directory

When running, Labelmerge automatically downloads and caches necessary conda envs to speed up subsequent runs.

By default, these are stored in the following directory:

```bash
~/.cache/labelmerge/
```

You can override this default cache location by setting the `LABELMERGE_CACHE_DIR` environment variable:

```bash
export LABELMERGE_CACHE_DIR=/path/to/custom/cache
```

This is useful when working on shared systems, when home directory storage is limited, or if you wish to isolate data per project or user.

## For Developers & Contributors

These steps are intended for people who want to contribute to the development of labelmerge or explore its internals.

### 1. Clone the labelmerge GitHub repository

```bash
git clone https://github.com/khanlab/labelmerge.git
cd labelmerge
```

---

### 2. Create and activate a new Conda environment

```bash
mamba env create -f labelmerge-dev.yml
eval "$(mamba shell hook --shell zsh)"
mamba activate labelmerge-dev
```

---

### 3. Run the development version of labelmerge

You can run labelmerge directly from the source directory using:

```bash
./labelmerge/run.py -h
```

This should print out the available command-line arguments for the tool.  
You’re now set up for development and contribution!

If you’re on an M-chip mac, prefix with CONDA_SUBDIR=osx-64 to ensure compatibility:

```bash
CONDA_SUBDIR=osx-64 ./labelmerge/run.py -h
```

---

## Troubleshooting

If you encounter issues while setting up labelmerge via Conda:

- Make sure you’re using the latest Conda:
  ```bash
  conda update -n base -c defaults conda
  ```
- Double-check that your environment is activated (`mamba activate labelmerge-env` or `labelmerge-dev`)
- Try creating a fresh environment if problems persist
- Search for similar issues or open a new one in the [GitHub issues](https://github.com/afids/autoafids/issues) page

---

Happy labelmerging!
