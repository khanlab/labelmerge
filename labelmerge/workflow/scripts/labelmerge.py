#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import xarray as xr


def load_atlas(atlas_path):
    """Loading relevant atlas data"""
    atlas = nib.load(atlas_path)
    data = atlas.get_fdata().astype(np.ushort)
    header = atlas.header
    affine = atlas.affine

    return data, header, affine


def assemble_mask(atlas, metadata, label, prefix=""):
    name = metadata[metadata["index"] == label].name.iloc[0]
    return (
        f"{prefix}{name}",
        xr.DataArray(atlas == label, dims=["x", "y", "z"]),
    )


def split_labels(atlas, metadata, prefix=""):
    return xr.Dataset(
        dict(
            [
                assemble_mask(atlas, metadata, label, prefix)
                for label in np.unique(atlas[atlas > 0])
            ]
        )
    )


def merge_labels(base_ds, overlay_ds):
    out_arr = np.zeros(list(dim for dim in base_ds.dims.values()))
    out_labels = list(range(len(base_ds) + len(overlay_ds)))
    out_names = out_labels.copy()

    overall_idx = 0
    for dataset in [base_ds, overlay_ds]:
        for name, arr in dataset.items():
            out_arr[arr] = overall_idx
            out_labels[overall_idx] = overall_idx + 1
            out_names[overall_idx] = name
            overall_idx += 1
    out_metadata = pd.DataFrame({"index": out_labels, "name": out_names})
    return out_arr, out_metadata


def gen_parser():
    parser = ArgumentParser()
    parser.add_argument("base_map")
    parser.add_argument("base_metadata")
    parser.add_argument("overlay_map")
    parser.add_argument("overlay_metadata")
    parser.add_argument("out_map_path")
    parser.add_argument("out_metadata_path")
    return parser


def main():
    parser = gen_parser()
    args = parser.parse_args()
    base_data, base_header, base_affine = load_atlas(Path(args.base_map))
    base_metadata = pd.read_csv(args.base_metadata, sep="\t")
    overlay_data, _, _ = load_atlas(Path(args.overlay_map))
    overlay_metadata = pd.read_csv(args.overlay_metadata, sep="\t")
    base_ds = split_labels(base_data, base_metadata, prefix="base ")
    overlay_ds = split_labels(overlay_data, overlay_metadata, prefix="overlay ")
    merged_map, merged_metadata = merge_labels(base_ds, overlay_ds)
    merged_img = nib.Nifti1Image(
        dataobj=merged_map, affine=base_affine, header=base_header
    )
    nib.save(merged_img, args.out_map_path)
    merged_metadata.to_csv(args.out_metadata_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
