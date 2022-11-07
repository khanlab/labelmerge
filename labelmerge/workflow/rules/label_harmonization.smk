from functools import partial
from pathlib import Path
from snakebids import bids


def aggregate_input(wildcards, checkpoint):
    checkpoint_output = checkpoint.get(**wildcards).output.binary_dir
    label_fname = Path(
        bids(label="{label_idx}", suffix="mask.nii.gz", **wildcards)
    ).name
    label_mask_path = str(Path(checkpoint_output) / label_fname)
    return expand(
        label_mask_path,
        label_idx=glob_wildcards(label_mask_path).label_idx,
    )





checkpoint split_base_labels:
    input:
        labelmap=inputs["labelmap"].input_path,
    params:
        bids_dir=config["bids_dir"],
    output:
        binary_dir=directory(
            str(
                Path(
                    bids(
                        root=str(Path(config["output_dir"]) / "labelmerge-work"),
                        suffix="dseg",
                        desc=config.get("base_desc", None),
                        **inputs["labelmap"].input_wildcards
                    )
                )
            )
        ),
    # TODO: Update container that has appropriate dependencies
    # container:
    #     "docker://khanlab/neuroglia-core"
    script:
        "../scripts/label_split.py"


checkpoint split_overlay_labels:
    input:
        labelmap=overlay["labelmap"].input_path,
    params:
        bids_dir=config["overlay_bids_dir"],
    output:
        binary_dir=directory(
            str(
                Path(
                    bids(
                        root=str(Path(config["output_dir"]) / "labelmerge-work"),
                        suffix="dseg",
                        desc=config.get("overlay_desc", None),
                        **overlay["labelmap"].input_wildcards
                    )
                )
            )
        ),
    # TODO: Update container that has appropriate dependencies
    # container:
    #     "docker://khanlab/neuroglia-core"
    script:
        "../scripts/label_split.py"


rule aggregate:
    input:
        base_labels=partial(aggregate_input, checkpoint=checkpoints.split_base_labels),
        overlay_labels=partial(aggregate_input, checkpoint=checkpoints.split_overlay_labels),
    output:
        touch(
            bids(
                root=str(Path(config["output_dir"]) / "aggregated"),
                suffix="aggregated",
                desc="combined",
                **inputs["labelmap"].input_wildcards
            )
        ),
    shell:
        "echo {input}"
