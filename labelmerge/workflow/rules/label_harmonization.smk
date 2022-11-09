from functools import partial
from pathlib import Path
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
        "base": expand(
            base_mask_path,
            label_idx=glob_wildcards(base_mask_path).label_idx,
        ),
        "overlay": expand(
            overlay_mask_path,
            label_idx=glob_wildcards(overlay_mask_path).label_idx,
        ),
    }



def choose_root(wildcards):
    return config["bids_dir"] if wildcards["desc"] == config["base_desc"] else config["overlay_bids_dir"]


def build_labelmap_path(wildcards):
    return {
        "labelmap": bids(
            root=choose_root(wildcards),
            datatype="anat",
            suffix="dseg.nii.gz",
            **inputs["labelmap"].input_wildcards
        ),
    }


checkpoint split_labels:
    input:
        unpack(build_labelmap_path)
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
