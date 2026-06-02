#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-tables", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--dataset-size", type=int, required=True)
    args = parser.parse_args()

    tables = []
    for table in args.full_tables:
        df = pd.read_csv(table, sep="\t")
        tables.append(df)

    data = pd.concat(tables, ignore_index=True)

    # Remove exact duplicated BUSCO hits
    dedup_columns = [
        "specie_taxid",
        "specie_name",
        "source",
        "data_id",
        "busco_id",
        "status",
        "sequence",
    ]

    before = len(data)

    data = data.drop_duplicates(
        subset=dedup_columns,
        keep="first"
    ).copy()

    after = len(data)

    print(f"Removed {before - after} duplicated BUSCO rows")

    required = {"specie_taxid", "specie_name", "busco_id"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    thresholds = range(40, 101, 10)
    summary_rows = []

    # Une espèce = taxid si disponible, sinon nom
    data["species_key"] = data["specie_taxid"].astype(str)
    data.loc[data["species_key"].isin(["", "nan", "None"]), "species_key"] = data["specie_name"]

    for busco_threshold in thresholds:
        # Complétude BUSCO par espèce
        species_busco_counts = (
            data.groupby("species_key")["busco_id"]
            .nunique()
            .reset_index(name="n_busco")
        )

        species_busco_counts["busco_complete_percentage"] = (
            species_busco_counts["n_busco"] / args.dataset_size * 100
        )

        selected_species = set(
            species_busco_counts.loc[
                species_busco_counts["busco_complete_percentage"] >= busco_threshold,
                "species_key",
            ]
        )

        data_species_filtered = data[data["species_key"].isin(selected_species)].copy()
        n_species_after_busco_filter = data_species_filtered["species_key"].nunique()

        for species_threshold in thresholds:
            if n_species_after_busco_filter == 0:
                filtered = data_species_filtered.iloc[0:0].copy()
                selected_buscos = set()
            else:
                busco_species_counts = (
                    data_species_filtered.groupby("busco_id")["species_key"]
                    .nunique()
                    .reset_index(name="n_species")
                )

                busco_species_counts["species_percentage"] = (
                    busco_species_counts["n_species"]
                    / n_species_after_busco_filter
                    * 100
                )

                selected_buscos = set(
                    busco_species_counts.loc[
                        busco_species_counts["species_percentage"] >= species_threshold,
                        "busco_id",
                    ]
                )

                filtered = data_species_filtered[
                    data_species_filtered["busco_id"].isin(selected_buscos)
                ].copy()

            out_file = (
                f"{args.output}_busco-c_{busco_threshold}"
                f"_gene-occupancy_{species_threshold}.tsv"
            )

            filtered.drop(columns=["species_key"], errors="ignore").to_csv(
                out_file,
                sep="\t",
                index=False,
            )

            summary_rows.append(
                {
                    "busco_complete_threshold": busco_threshold,
                    "species_presence_threshold": species_threshold,
                    "n_species": filtered["species_key"].nunique()
                    if not filtered.empty
                    else 0,
                    "n_orthogroups": filtered["busco_id"].nunique()
                    if not filtered.empty
                    else 0,
                    "n_rows": len(filtered),
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary_file = f"{args.output}_summary.tsv"
    summary.to_csv(summary_file, sep="\t", index=False)

        # Plot 1: heatmap orthogroups with species annotation
    pivot_og = summary.pivot(
        index="species_presence_threshold",
        columns="busco_complete_threshold",
        values="n_orthogroups",
    )

    pivot_sp = summary.pivot(
        index="species_presence_threshold",
        columns="busco_complete_threshold",
        values="n_species",
    )

    fig, ax = plt.subplots(figsize=(9, 7))

    im = ax.imshow(pivot_og.values, aspect="auto", origin="lower")

    ax.set_xticks(range(len(pivot_og.columns)))
    ax.set_xticklabels(pivot_og.columns)
    ax.set_yticks(range(len(pivot_og.index)))
    ax.set_yticklabels(pivot_og.index)

    ax.set_xlabel("BUSCO completeness threshold per species (%)")
    ax.set_ylabel("BUSCO presence threshold across species (%)")
    ax.set_title("Orthogroups retained across threshold combinations")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Number of orthogroups")

    for i in range(pivot_og.shape[0]):
        for j in range(pivot_og.shape[1]):
            og = int(pivot_og.iloc[i, j])
            sp = int(pivot_sp.iloc[i, j])
            ax.text(j, i, f"{og}\n({sp} sp.)", ha="center", va="center", fontsize=7)

    plt.tight_layout()

    for ext in ["pdf", "png", "svg"]:
        plt.savefig(f"{args.output}_heatmap.{ext}", dpi=300)

    plt.close()

    # Plot 2: compromise scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))

    scatter = ax.scatter(
        summary["n_species"],
        summary["n_orthogroups"],
        c=summary["busco_complete_threshold"],
        s=summary["species_presence_threshold"] * 2,
        alpha=0.8,
    )

    ax.set_xlabel("Number of species retained")
    ax.set_ylabel("Number of orthogroups retained")
    ax.set_title("Species / orthogroups trade-off")

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("BUSCO completeness threshold (%)")

    plt.tight_layout()

    for ext in ["pdf", "png", "svg"]:
        plt.savefig(f"{args.output}_tradeoff.{ext}", dpi=300)

    plt.close()


if __name__ == "__main__":
    main()