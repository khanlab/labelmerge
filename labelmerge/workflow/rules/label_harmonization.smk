from pathlib import Path
from snakebids import bids


checkpoint split_labels:
    input:
        labelmap=inputs["labelmap"].input_path,
    output:
        binary_dir=directory(
            str(
                Path(
                    bids(
                        root=str(Path(config["output_dir"]) / "labelmerge-work"),
                        suffix="dseg",
                        **inputs["labelmap"].input_wildcards
                    )
                )
            )
        ),
    # Update container that has appropriate dependencies
    # container:
    #     "docker://khanlab/neuroglia-core"
    script:
        '../scripts/label_split.py'

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_labels.get(**wildcards).output.binary_dir
    label_fname = Path(
            bids(label="{label_idx}", suffix="dseg.nii.gz", **wildcards)
        ).name
    label_mask_path = str(Path(checkpoint_output) / label_fname)
    return expand(
        label_mask_path,
        label_idx=glob_wildcards(label_mask_path).label_idx,
    )


rule aggregate:
    input:
        aggregate_input,
    output:
        touch(
            bids(
                root=str(Path(config["output_dir"]) / "aggregated"),
                suffix="aggregated",
                **inputs["labelmap"].input_wildcards
            )
        ),
    shell:
        "echo {input}"
