#!/usr/bin/env python3

import argparse
import pandas as pd


SOURCE_PRIORITY = {
    "user_proteomes": 1,
    "uniprot_proteomes": 2,
    "ncbi_proteomes": 3,
    "ncbi_assemblies": 4,
    "tsa_transcriptomes": 5,
    "sra_transcriptomes": 5,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-tables", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    dfs = []

    for table in args.full_tables:
        df = pd.read_csv(table, sep="\t")
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    df["source_rank"] = df["source"].map(SOURCE_PRIORITY).fillna(999)

    df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(0)
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(0)

    df = df.sort_values(
        by=["busco_id", "source_rank", "length", "score"],
        ascending=[True, True, False, False],
    )

    best = df.groupby("busco_id", as_index=False).first()

    best = best.drop(columns=["source_rank"])

    out_fields = [
        "specie_taxid",
        "specie_name",
        "source",
        "data_id",
        "busco_id",
        "status",
        "sequence",
        "score",
        "length",
        "orthodb_url",
        "description",
    ]

    best = best[out_fields]

    best.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()