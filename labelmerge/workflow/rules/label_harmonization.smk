from functools import partial
import glob
from pathlib import Path
from bids.layout import parse_file_entities
from snakebids import bids


def split_labels_path(wildcards):
    return str(
        Path(checkpoints.split_labels.get(**wildcards).output.binary_dir)
        / Path(bids(label="{label_idx}", suffix="mask.nii.gz", **wildcards)).name
    )


def aggregate_inputs(wildcards):
    base_mask_path = split_labels_path(dict(wildcards, desc=config["base_desc"]))
    overlay_mask_path = split_labels_path(dict(wildcards, desc=config["overlay_desc"]))
    return {
        "base_masks": expand(
            base_mask_path,
            label_idx=glob_wildcards(base_mask_path).label_idx,
        ),
        "base_metadata": build_metadata_path(
            dict(wildcards, desc=config["base_desc"])
        ),
        "overlay_masks": expand(
            overlay_mask_path,
            label_idx=glob_wildcards(overlay_mask_path).label_idx,
        ),
        "overlay_metadata": build_metadata_path(
            dict(wildcards, desc=config["overlay_desc"])
        ),
    }


def choose_root(wildcards):
    return (
        config["bids_dir"]
        if wildcards["desc"] == config["base_desc"]
        else config["overlay_bids_dir"]
    )


def build_labelmap_path(wildcards):
    return bids(
        root=choose_root(wildcards),
        datatype="anat",
        suffix="dseg.nii.gz",
        **inputs["labelmap"].input_wildcards,
    )


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

    return tsv_file


def build_metadata_path(wildcards):
    atlas_path = expand(build_labelmap_path(wildcards), **wildcards)[0]
    return load_metadata(Path(atlas_path), choose_root(wildcards))


checkpoint split_labels:
    input:
        labelmap=build_labelmap_path,
        metadata=build_metadata_path,
    params:
        bids_dir=choose_root,
    output:
        binary_dir=directory(
            bids(
                root=str(Path(config["output_dir"]) / "labelmerge-work"),
                suffix="dseg",
                **inputs["labelmap"].input_wildcards
            )
        ),
    # TODO: Update container that has appropriate dependencies
    # container:
    #     "docker://khanlab/neuroglia-core"
    script:
        "../scripts/label_split.py"


rule aggregate:
    input:
        unpack(aggregate_inputs),
    output:
        touch(
            bids(
                root=str(Path(config["output_dir"]) / "aggregated"),
                suffix="aggregated",
                desc="combined",
                **dict(
                    item
                    for item in inputs["labelmap"].input_wildcards.items()
                    if item[0] != "desc"
                )
            )
        ),
    shell:
        "echo {input}"
