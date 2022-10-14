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
    container:
        "docker://khanlab/neuroglia-core"
    shell:
        "mkdir -p {output.binary_dir}; c3d {input} -split -oo {output.binary_dir}/label%02d.nii.gz"


def aggregate_input(wildcards):
    checkpoint_output = checkpoints.split_labels.get(**wildcards).output.binary_dir
    label_mask_path = str(Path(checkpoint_output) / "label{i}.nii.gz")
    return expand(
        label_mask_path,
        i=glob_wildcards(label_mask_path).i,
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
