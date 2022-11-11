import re

import pandas as pd


ALPHANUMERIC_RE = re.compile(r"[^a-zA-Z\d]")


def read_table(table_path):
    label_df = pd.read_table(table_path)
    label_df = label_df[
        ["index", "name"] + (["abbreviation"] if "abbreviation" in label_df else [])
    ]
    print(label_df)
    if "abbreviation" not in label_df:
        label_df["abbreviation"] = label_df["name"]
    label_df = label_df.assign(
        abbreviation=label_df["abbreviation"].str.replace(
            ALPHANUMERIC_RE, "", regex=True
        )
    )
    return label_df


def assemble_config_tsv(base_path, overlay_path):
    return pd.concat(
        [
            read_table(base_path).assign(priority=1).assign(source="base"),
            read_table(overlay_path).assign(priority=2).assign(source="overlay"),
        ]
    )


def main():
    config_df = assemble_config_tsv(
        snakemake.input["base_metadata"], snakemake.input["overlay_metadata"]
    )
    config_df.to_csv(snakemake.output["config_tsv"], sep="\t", index=False)


if __name__ == "__main__":
    main()
