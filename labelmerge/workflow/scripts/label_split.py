#!/usr/bin/env python
import glob
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


def load_metadata(atlas_path, bids_dir, smk_wildcards):
    """Load metadata associated with atlas"""
    # Walk backwards to find dseg until bids_dir is hit
    tsv_file = None
    cur_dir = atlas_path.parent
    while cur_dir != Path(bids_dir).parent and not tsv_file:
        # Check number of dseg files found
        dseg_files = list(glob.iglob(f"{cur_dir}/*dseg.tsv"))
        num_dsegs = len(dseg_files)

        if num_dsegs == 1:
            # If single file, assume association
            tsv_file = Path(dseg_files[0])
        elif num_dsegs > 1:
            # If multiple files, assume desc entity exists
            for dseg_file in dseg_files:
                if smk_wildcards["desc"] in dseg_file:
                    tsv_file = dseg_file
        else:
            # Move up a directory
            cur_dir = cur_dir.parent

    # If still no file found
    if not tsv_file:
        raise FileNotFoundError("No associated tsv file found")

    # Read associated file and create new column storing description
    metadata = pd.read_csv(tsv_file, sep="\t")
    metadata["BIDS Name"] = metadata["Name"].str.title().str.replace(" ", "")

    return metadata


def label_split(atlas_path, output_dir, smk_wildcards, bids_dir):
    """Split labels from atlas into individual files"""

    # Create parent directory
    Path(output_dir).mkdir(parents=True)

    # Load atlas + metadata
    atlas_path = Path(atlas_path)
    atlas_data, atlas_header, atlas_affine = load_atlas(atlas_path)
    atlas_metadata = load_metadata(atlas_path, bids_dir, smk_wildcards)

    # Extract & save unique labels
    for label in np.unique(atlas_data[atlas_data > 0]):
        # Get indices for label and create a binary mask
        label_indices = np.argwhere(atlas_data == label)
        label_mask = np.ma.make_mask(atlas_data[label_indices], copy=True)
        label_img = nib.Nifti1Image(
            dataobj=label_mask, affine=atlas_affine, header=atlas_header
        )

        # Grab label description
        try:
            label_desc = atlas_metadata[atlas_metadata["Index"] == int(label)][
                "BIDS Name"
            ].values[0]
        except IndexError:
            raise ValueError(
                f"Label {int(label)} does not exist in the associated tsv file"
            )

        # Set file name
        label_fname = Path(
            bids(label=label_desc, suffix="mask.nii.gz", **smk_wildcards)
        ).name

        nib.save(label_img, f"{output_dir}/{label_fname}")


if __name__ == "__main__":
    label_split(
        atlas_path=snakemake.input["labelmap"],  # noqa: F821
        output_dir=snakemake.output["binary_dir"],  # noqa: F821
        smk_wildcards=snakemake.wildcards,  # noqa: F821
        bids_dir=snakemake.params["bids_dir"],  # noqa: F821
    )
