# Contributing to Labelmerge

Labelmerge python package dependencies are managed with Poetry (`v1.2.0+`), which
you will need installed on your machine. You can find installation instructions 
on the [Poetry website](https://python-poetry.org/docs/master/#installation).

_Note: These instructions are only recommended if you are making changes to the
Labelmerge codebase and committing these back to the repository, or if you are 
using Snakemake's cluster execution profiles. If not, it is easier to run 
Labelmerge using the packaged singularity container (e.g. 
`docker://khanlab/labelmerge:latest`)._

## Setup the development environment

Clone the repository and install all dependencies (including `dev`) with poetry:

```
git clone https://github.com/khanlab/labelmerge.git 
cd labelmerge
poetry install --with dev 
```

Poetry will automatically create a virtual environment. To customize where 
these virtual environments are stored, see the poetry docs 
[here](https://python-poetry.org/docs/configuration/)

Then, you can run labelmerge with:

```
poetry run labelmerge
```

or you can activate a virtual environment shell and run labelmerge directly:

```
poetry shell
labelmerge
```

You can exit the poetry shell with `exit`

## Running and fixing code format quality

labelmerge uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner.
You can see what commands are available by running:

```
poetry run poe 
```

We use a a few tools, including `black`, `flake8`, `isort`, `snakefmt`, and 
`yamlfix` to ensure formatting and style of our codebase is consistent. There 
are two task runners you can use to check and fix your code, which can be 
invoked with:

```
poetry run poe quality-check
poetry run poe quality
```

_Note: If you are in a poetry shell, you do not need to prepend `poetry run` to
the command._

## Dry-run / testing your workflow
Using Snakemake\'s dry-run option (`--dry-run`/`-n`) is an easy way to verify
any changes made to the workflow are working direcctly. The `test/data` folder 
contains a _fake* BIDS dataset (i.e. dataset with zero-sized files) that are 
useful for verifying different aspects of the workflow. These dry-run tests are 
part of the automated Github actions that are run for every commit.

You can invoke the pre-configured task via 
[poethepoet](https://github.com/nat-n/poethepoet) to perform a dry-run:

```
poetry run poe test
```

This performs a number of tests, involving different scenarios in which a user
may use labelmerge.