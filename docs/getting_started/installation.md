# Installation

BIDS App for combining multiple atlas parcellations

## Requirements

* Docker (Intel Mac/Windows/Linux) or Singularity (Linux)

### Notes
* Inputs to Labelmerge should be organized as [BIDS derivatives](https://bids-specification.readthedocs.io/en/stable/01-introduction.html), taking one "base" and one "overlay" atlas

## Docker on Windows / Mac (Intel) / Linux

The Labelmerge BIDS App is available on DockerHub as versioned releases.
Instructions can be found in the [Docker](https://labelmerge.readthedocs.io/en/stable/getting_started/docker.html) documentation page.

### Pros
* Compatible with non-Linux systems
* All dependencies are in a single container

### Cons
* Typically not possible on shared machines
* Cannot use Snakemake cluster execution profiles
* Cannot edit code

## Singularity Container

The same Docker container can also be used with Singularity (now Apptainer).
Instructions can be found in the [Singularity / Apptainer](https://labelmerge.readthedocs.io/en/stable/getting_started/singularity.html) documentation page.

### Pros
* All dependencies are in a single container, stored as a single file (.sif)
* Compatible on shared systems with Singularity installed

### Cons
* Cannot use Snakemake cluster execution profiles
* Cannot edit code

## Python Environment with Singularity Dependencies

Instructions can be found in the [Contributing](https://labelmerge.readthedocs.io/en/stable/contributing/contributing.html) documentation page.

### Pros
* Flexibility to modify code

### Cons
* Only compatible on systems with Singularity for external dependencies

## Pip Install
`labelmerge` can also be installed using pip:

```bash
pip install labelmerge
```

