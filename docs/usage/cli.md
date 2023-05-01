# Command Line Interface (CLI)

## labelmerge CLI
The following can also be seen by entering `labelmerge -h` into your terminal.

These are all the required and optional arguments labelmerge accepts in order to 
run flexibly on many different input data types and with many options. In most 
cases, only the required arguments are needed.

```{argparse}
---
filename: ../labelmerge/run.py
func: get_parser
prog: labelmerge
---
```

## Snakemake CLI
In addition to the above command line arguments, Snakemake arguments can also be
passed at the labelmerge command line.

The complete list of [Snakemake](https://snakemake.readthedocs.io/en/stable/) 
arguments are below, and most act to determine your environment and app
behaviours. They will likely only need to be used for running in cloud
environments or troubleshooting. These can be listed from the command line with
`scatter --help-snakemake`.

```{argparse}
---
module: snakemake
func: get_argument_parser
prog: snakemake
---
```
