## Labelmerge
![Version](https://img.shields.io/github/v/tag/khanlab/labelmerge?label=version)
![Python3](https://img.shields.io/badge/python-3.8_|_3.9_|_3.10-blue.svg)
![Docker Pulls](https://img.shields.io/docker/pulls/khanlab/labelmerge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7636410.svg)](https://doi.org/10.5281/zenodo.7636410)

Labelmerge: A BIDS app for merging multiple label maps and metadata into a single parcellation scheme.

Examination of the brain has led to the creation and distribution of numerous parcellation schemes, studying various features or aspects including gross anatomy, cytoarchitecture, myeloarchitecture, functional connectivty, or structural connectivity. With the growing number of available schemes, there exists a need for methods to robustly merge different these different atlases from numerous sources depending on the aims of the specific study.  

Currently, combining multiple atlases requires image-processing tools to manually remove and add the regions of interest from the respective atlases. While effective, this task becomes inefficient and is prone to error when looking to combine specific labels of interest from multiple atlases or across a group of subjects. To that end we developed Labelmerge, a Brain Imaging Data Structure (BIDS) app that combines the parcellations of multiple atlases into a single parcellation scheme that can be applied towards downstream analysis.


### Contributing
Clone the git repository. Labelmerge dependencies are managed with Poetry
(version 1.2.x), which you'll need installed on your machine.
You can find instructions on the Poetry
[website](https://python-poetry.org/docs/).

Then, setup the development environment with the following commands:

```
poetry install
poetry run poe setup
```

Labelmerge uses poethepoet as a task runner.
You can see what commands are available by running:

```
poetry run poe
```

If you wish, you can also run poe [command] directly by installing poethepoet
on your system. Follow the install instructions at the link above.

Labelmerge uses pre-commit hooks (installed via the poe setup command above) to 
lint and format code (we use black, isort, flake8). By default, these hooks are
run on every commit.

Please be sure they all pass before making a PR.


