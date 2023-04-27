from functools import partial
import glob
from pathlib import Path
import re

from bids.layout import parse_file_entities
from snakebids import bids, generate_inputs


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
    base_metadata = load_metadata(
        Path(expand(base_inputs["labelmap"].path, **wildcards)[0]),
        config["bids_dir"],
    )
    overlay_metadata = load_metadata(
        Path(expand(overlay_inputs["labelmap"].path, **wildcards)[0]),
        config["overlay_bids_dir"],
    )
    return {
        "base_metadata": base_metadata,
        "overlay_metadata": overlay_metadata,
    }


rule merge_labels:
    input:
        unpack(build_metadata_path),
        base_map=base_inputs["labelmap"].path,
        overlay_map=overlay_inputs["labelmap"].path,
    output:
        merged_map=bids(
            root=str(Path(config["output_dir"]) / "combined"),
            suffix="dseg.nii.gz",
            desc="combined",
            **base_inputs["labelmap"].wildcards,
        ),
        merged_metadata=bids(
            root=str(Path(config["output_dir"]) / "combined"),
            suffix="dseg.tsv",
            desc="combined",
            **base_inputs["labelmap"].wildcards,
        ),
    params:
        base_exceptions=f"--base_exceptions {' '.join(config['base_exceptions'])}"
        if config.get("base_exceptions")
        else "",
        overlay_exceptions=f"--overlay_exceptions {' '.join(config['overlay_exceptions'])}"
        if config.get("overlay_exceptions")
        else "",
        base_drops=f"--base_drops {' '.join(config['base_drops'])}"
        if config.get("base_drops")
        else "",
        overlay_drops=f"--overlay_drops {' '.join(config['overlay_drops'])}"
        if config.get("overlay_drops")
        else "",
    resources:
        script=str(Path(workflow.basedir) / "scripts" / "labelmerge.py"),
    shell:
        "python3 {resources.script} {input.base_map} {input.base_metadata} "
        "{input.overlay_map} {input.overlay_metadata} "
        "{output.merged_map} {output.merged_metadata} "
        "{params.base_exceptions} {params.overlay_exceptions} "
        "{params.base_drops} {params.overlay_drops}"
