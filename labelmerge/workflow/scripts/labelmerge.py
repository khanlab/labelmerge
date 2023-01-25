#!/usr/bin/env python
"""Script to merge two labelmaps with their BIDS metadata."""
from __future__ import annotations

from argparse import ArgumentParser
from os import PathLike
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import xarray as xr
from numpy.typing import ArrayLike


def load_atlas(atlas_path: PathLike):
    """Loading relevant atlas data"""
    atlas = nib.load(atlas_path)
    data = atlas.get_fdata().astype(np.int_)
    header = atlas.header
    affine = atlas.affine

    return data, header, affine


class MetadataError(Exception):
    """Error raised when there's an issue with the metadata."""


def assemble_mask(
    atlas: ArrayLike,
    metadata: pd.DataFrame,
    label: int,
    prefix: str = "",
):
    try:
        name: str = metadata[metadata["index"] == label].name.iloc[0]
    except IndexError as err:
        raise MetadataError(
            f"Label with index {label} from {prefix}atlas not found in "
            "metadata table."
        ) from err
    return (
        f"{prefix}{name}",
        xr.DataArray(atlas == label, dims=["x", "y", "z"]),
    )


def split_labels(
    atlas: np.ndarray,
    metadata: pd.DataFrame,
    prefix: str = "",
    exceptions: list[int] | None = None,
    drops: list[int] | None = None,
) -> list[xr.Dataset]:
    if exceptions is None:
        exceptions = []
    if drops is None:
        drops = []
    unique_vals = np.unique(atlas[atlas > 0])
    all_labels: pd.Series[int] = metadata["index"]
    if not set(unique_vals) <= set(all_labels):
        unlabeled_vals = ", ".join(
            str(val) for val in set(unique_vals) - set(all_labels)
        )
        raise MetadataError(
            f"Labels with indices {unlabeled_vals} from {prefix}atlas not "
            "found in metadata table"
        )
    normal_ds = xr.Dataset(
        dict(
            [
                assemble_mask(atlas, metadata, label, prefix)
                for label in unique_vals
                if label not in exceptions + drops
            ]
        )
    )
    if not exceptions:
        return [normal_ds]
    exception_ds = xr.Dataset(
        dict(
            [
                assemble_mask(atlas, metadata, label, prefix)
                for label in all_labels
                if label in exceptions
            ]
        )
    )
    return [normal_ds, exception_ds]


def merge_labels(datasets: list[xr.Dataset]):
    """Merge several labelmaps into one output labelmap with metadata.

    Parameters
    ----------
    datasets: list of Dataset
        A list of datasets where each variable name is a label name and each
        DataArray is the corresponding mask. Earlier entries in the list of
        datasets are overwritten by later entries.

    Returns
    -------
    ndarray
        A 3D labelmap with the merged label data.
    DataFrame
        The metadata table containing label names.
    """
    out_arr: np.ndarray = np.zeros(
        list(dim for dim in datasets[0].dims.values()), dtype=np.int_
    )
    out_labels = [
        idx + 1 for idx in range(sum(len(dataset) for dataset in datasets))
    ]
    out_names = ["" for _ in out_labels]

    overall_idx = 0
    for dataset in datasets:
        for name, arr in dataset.items():
            out_arr[arr] = int(overall_idx + 1)
            out_names[overall_idx] = name
            overall_idx += 1
    out_metadata = pd.DataFrame({"index": out_labels, "name": out_names})
    return out_arr, out_metadata


def gen_parser() -> ArgumentParser:
    """Generate the CLI parser for this script."""
    parser = ArgumentParser()
    parser.add_argument("base_map")
    parser.add_argument("base_metadata")
    parser.add_argument("overlay_map")
    parser.add_argument("overlay_metadata")
    parser.add_argument("out_map_path")
    parser.add_argument("out_metadata_path")
    parser.add_argument(
        "--base_exceptions",
        nargs="*",
        help=(
            "Space separated list of integer labels from the base image to "
            "keep over overlay labels at the same voxels."
        ),
        type=int,
    )
    parser.add_argument(
        "--overlay_exceptions",
        nargs="*",
        help=(
            "Space separated list of integer labels from the overlay image to "
            "be overwritten by labels from the base image."
        ),
        type=int,
    )
    parser.add_argument(
        "--base_drops",
        nargs="*",
        help=(
            "Space separated list of integer labels from the base image to "
            "drop from the output map."
        ),
        type=int,
    )
    parser.add_argument(
        "--overlay_drops",
        nargs="*",
        help=(
            "Space separated list of integer labels from the overlay image to "
            "drop from the output map."
        ),
        type=int,
    )
    return parser


def main():
    parser = gen_parser()
    args = parser.parse_args()
    base_data, base_header, base_affine = load_atlas(Path(args.base_map))
    base_metadata = pd.read_csv(args.base_metadata, sep="\t")
    overlay_data, _, _ = load_atlas(Path(args.overlay_map))
    overlay_metadata = pd.read_csv(args.overlay_metadata, sep="\t")
    base_exceptions = args.base_exceptions if args.base_exceptions else []
    overlay_exceptions = (
        args.overlay_exceptions if args.overlay_exceptions else []
    )
    base_drops = args.base_drops if args.base_drops else []
    overlay_drops = args.overlay_drops if args.overlay_drops else []
    base_datasets = split_labels(
        base_data,
        base_metadata,
        prefix="base ",
        exceptions=base_exceptions,
        drops=base_drops,
    )
    overlay_datasets = split_labels(
        overlay_data,
        overlay_metadata,
        prefix="overlay ",
        exceptions=overlay_exceptions,
        drops=overlay_drops,
    )
    # Note that overlay exceptions are ignored
    merged_map, merged_metadata = merge_labels(
        overlay_datasets[1:]
        + base_datasets[0:1]
        + overlay_datasets[0:1]
        + base_datasets[1:]
    )
    merged_img = nib.Nifti1Image(
        dataobj=merged_map, affine=base_affine, header=base_header
    )
    nib.save(merged_img, args.out_map_path)
    merged_metadata.to_csv(args.out_metadata_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
