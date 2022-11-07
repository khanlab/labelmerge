#!/usr/bin/env python
import glob
from pathlib import Path
from sys import exit

from bids.layout import parse_file_entities
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


def load_metadata(atlas_path, bids_dir):
    """Load metadata associated with atlas"""
    # Walk backwards to find dseg until bids_dir is hit
    tsv_file = None
    cur_dir = atlas_path.parent
    atlas_entities = parse_file_entities(atlas_path)
    del atlas_entities["extension"]
    while cur_dir != Path(bids_dir).parent and tsv_file is None:
        # Check number of dseg files found
        dseg_files = list(glob.iglob(f"{cur_dir}/*dseg.tsv"))

        for candidate in dseg_files:
            possible = True
            candidate_entities = parse_file_entities(candidate)
            del candidate_entities["extension"]
            for key, val in candidate_entities.items():
                if (key not in atlas_entities) or (atlas_entities[key] != val):
                    possible = False
            if possible:
                if tsv_file is None:
                    tsv_file = candidate
                else:
                    raise ValueError(
                        "Multiple applicable metadata files applicable at this level "
                        "of the directory hierarchy"
                    )
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
    atlas_metadata = load_metadata(atlas_path, bids_dir)

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
        except:
            raise ValueError(
                f"Label {int(label)} does not exist in the associated tsv file"
            )
            exit(1)

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
