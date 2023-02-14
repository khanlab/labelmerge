## Labelmerge
![Version](https://img.shields.io/github/v/tag/khanlab/labelmerge?label=version)
![Python3](https://img.shields.io/badge/python-3.8_|_3.9_|_3.10-blue.svg)
![Docker Pulls](https://img.shields.io/docker/pulls/khanlab/labelmerge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7636410.svg)](https://doi.org/10.5281/zenodo.7636410)

Merge multiple label maps.

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


