# Installation

## Requirements

Labelmerge can be run either with Conda, or with containers
-   Conda (Linux/Mac only, no containers needed)
-   **--OR--**  Docker (Intel Mac/Windows/Linux) or Apptainer (formerly known as Singularity) (Linux)
-   For those wishing to contribute or modify the code, see [Contributing to Labelmerge](https://labelmerge.readthedocs.io/en/latest/contributing/contributing.html).
-   GPU not required

### Notes
* Inputs to Labelmerge should be organized as [BIDS derivatives](https://bids-specification.readthedocs.io/en/stable/05-derivatives/03-imaging.html#segmentations), taking one "base" and one "overlay" atlas


## Comparison of methods for running AutoAFIDs

There are several different ways of running AutoAFIDs:

1. Conda Environment (Linux/macOS)
2. Apptainer/Singularity Container on Linux
3. Docker Container on Windows/Mac (Intel)/Linux

## Conda

Labelmerge is available via [Conda](https://docs.conda.io/), offering a container-free way to run the tool on Linux and macOS systems.

### Pros:
- No need for any containers
- Easy to install via `conda`
- Compatible with Snakemake execution profiles
- Good option for users who prefer Python virtual environments

### Cons:
- Not compatible with Windows
- Requires installation of Conda or Miniconda

## Docker on Windows/macOS/Linux

The Labelmerge BIDS App is available as a Docker container. Instructions can be found in the [Docker](https://labelmerge.readthedocs.io/en/latest/getting_started/docker.html) documentation page.

### Pros
- Compatible with non-Linux systems, including Apple Silicon (M1, M2)
- All dependencies are contained within a single container

### Cons
- Typically not possible on shared machines
- Cannot use Snakemake cluster execution profiles
- Cannot edit code within the container

## Singularity Container

The same Docker container can also be used with Singularity (now Apptainer). Instructions can be found in the [Singularity](https://labelmerge.readthedocs.io/en/latest/getting_started/singularity.html) documentation page.

### Pros
- All dependencies are packaged within a single `.sif` file
- Compatible with shared systems where Singularity is installed

### Cons
- Cannot use Snakemake cluster execution profiles
- Cannot edit code within the container

