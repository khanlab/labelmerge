#!/usr/bin/env python
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from snakebids import bids


def load_atlas(atlas_path):
    """Loading relevant atlas data"""
    atlas = nib.load(atlas_path)
    data = atlas.get_fdata().astype(np.ushort)
    header = atlas.header
    affine = atlas.affine

    return data, header, affine


def load_metadata(atlas_path):
    """Load metadata associated with atlas"""
    # Find and read associated metadata file with atlas
    try:
        tsv_path = atlas_path.replace("nii.gz", "tsv")
        metadata = pd.read_csv(tsv_path, sep="\t")
    except:
        raise FileNotFoundError(f"{tsv_path} file not found")

    # Create new column storing description for filename
    metadata["BIDS Name"] = metadata["Name"].str.title().str.replace(" ", "")

    return metadata


def label_split(atlas_path, output_dir, smk_wildcards):
    """Split labels from atlas into individual files"""

    # Create parent directory
    Path(output_dir).mkdir(parents=True)

    # Load atlas + metadata
    atlas_data, atlas_header, atlas_affine = load_atlas(atlas_path)
    atlas_metadata = load_metadata(atlas_path)

    # Extract & save unique labels
    for label in np.unique(atlas_data[atlas_data > 0]):
        # Get indices for label and create a binary mask
        label_indices = np.argwhere(atlas_data == label)
        label_mask = np.ma.make_mask(atlas_data[label_indices], copy=True)
        label_img = nib.Nifti1Image(
            dataobj=label_mask, affine=atlas_affine, header=atlas_header
        )

        # Grab label description
        label_desc = atlas_metadata[atlas_metadata["Index"] == int(label)][
            "BIDS Name"
        ].values[0]

        # Set file name
        label_fname = Path(
            bids(label=label_desc, suffix="mask.nii.gz", **smk_wildcards)
        ).name

        nib.save(label_img, f"{output_dir}/{label_fname}")


if __name__ == "__main__":
    label_split(
        atlas_path=snakemake.input["labelmap"],
        output_dir=snakemake.output["binary_dir"],
        smk_wildcards=snakemake.wildcards,
    )
