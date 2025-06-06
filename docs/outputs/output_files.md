In the specified `/path/to/output/dir`, there will be the following outputs:

```
/path/to/output/dir/
└── config
├── combined
├── .snakebids
└── .snakemake
```

The `config` folder, along with the hidden `.snakebids` and `.snakemake` folders
contain a record of the code and parameters used, and paths to the inputs.

## combined Directory 
After running the workflow, the `/path/to/output/dir` folder will contain a `combined` directory. The combined label nii.gz image and the combined .tsv file will be in the `combined` directory with the following structure:

```
combined/
└── *_desc-combined_*.nii.gz
└── *_desc-combined_*.tsv
```
