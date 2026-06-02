#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SOURCES = {
    "proteomes_uniprot": "fetch_uniprot_proteomes/*.uniprot_proteome_ids.tsv",
    "ncbi_assemblies": "fetch_ncbi_assemblies/*.ncbi_accessions.tsv",
    "ncbi_proteomes": "fetch_ncbi_assemblies/*.ncbi_accessions_for_proteins.tsv",
    "tsa_transcriptomes": "fetch_tsa_transcriptomes/*.tsa_ids.tsv",
    "sra_runs": "fetch_sra_reads/*.sra_runs.tsv",
}


def clean_species_name(name: str) -> str:
    return (
        name.strip()
        .lower()
        .replace(" ", "-")
        .replace(".", "-")
    )


def read_ids(cache_dir: Path):
    data = defaultdict(lambda: defaultdict(set))
    species_names = {}

    for source, pattern in SOURCES.items():
        for file in cache_dir.glob(pattern):
            with open(file, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")

                for row in reader:
                    taxid = row.get("taxid", "").strip()
                    specie = clean_species_name(row.get("specie", ""))
                    data_id = None

                    for col in ["sra", "proteome", "genome", "tsa"]:
                        if col in row and row[col].strip():
                            data_id = row[col].strip()
                            break

                    if not taxid or not specie or not data_id:
                        continue

                    key = (taxid, specie)
                    species_names[key] = specie
                    data[key][source].add(data_id)

    return data, species_names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache-dir",
        default="/shared/projects/metainvert/phyloconstructor2/cache",
    )
    parser.add_argument("--output-prefix", default="available_data")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    data, species_names = read_ids(cache_dir)

    columns = [
        "taxid",
        "specie",
        "proteomes_uniprot",
        "ncbi_assemblies",
        "ncbi_proteomes",
        "tsa_transcriptomes",
        "sra_runs",
    ]

    ids_file = f"{args.output_prefix}_ids.tsv"
    counts_file = f"{args.output_prefix}_counts.tsv"

    with open(ids_file, "w", newline="") as out_ids, open(counts_file, "w", newline="") as out_counts:
        ids_writer = csv.DictWriter(out_ids, fieldnames=columns, delimiter="\t")
        counts_writer = csv.DictWriter(out_counts, fieldnames=columns, delimiter="\t")

        ids_writer.writeheader()
        counts_writer.writeheader()

        for taxid, specie in sorted(data.keys(), key=lambda x: (x[1], x[0])):
            ids_row = {
                "taxid": taxid,
                "specie": specie,
            }
            counts_row = {
                "taxid": taxid,
                "specie": specie,
            }

            for source in columns[2:]:
                values = sorted(data[(taxid, specie)].get(source, []))
                ids_row[source] = ",".join(values)
                counts_row[source] = len(values)

            ids_writer.writerow(ids_row)
            counts_writer.writerow(counts_row)

    print(f"Wrote {ids_file}")
    print(f"Wrote {counts_file}")


if __name__ == "__main__":
    main()