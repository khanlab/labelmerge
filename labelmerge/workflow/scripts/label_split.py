#!/usr/bin/env python
from pathlib import Path

import nibabel as nib
import numpy as np
from snakebids import bids


def load_atlas(atlas_path):
    """Loading relevant atlas data"""
    atlas = nib.load(atlas_path)
    data = atlas.get_fdata()
    header = atlas.header
    affine = atlas.affine

    return data, header, affine


def label_split(atlas_path, output_dir, smk_wildcards):
    """Split labels from atlas into individual files"""

    # Load atlas
    atlas_data, atlas_header, atlas_affine = load_atlas(atlas_path)
    Path(output_dir).mkdir(parents=True)

    # Extract & save unique labels
    for label in np.unique(atlas_data[atlas_data > 0]):
        # Get indices for label and create a binary mask
        label_indices = np.argwhere(atlas_data == label)
        label_mask = np.ma.make_mask(atlas_data[label_indices], copy=True)

        label_img = nib.Nifti1Image(
            dataobj=label_mask, affine=atlas_affine, header=atlas_header
        )

        # Set file name
        label_fname = Path(
            bids(label=int(label), suffix="dseg.nii.gz", **smk_wildcards)
        ).name

        nib.save(label_img, f"{output_dir}/{label_fname}")


if __name__ == "__main__":
    label_split(
        atlas_path=snakemake.input["labelmap"],
        output_dir=snakemake.output["binary_dir"],
        smk_wildcards=snakemake.wildcards,
    )
