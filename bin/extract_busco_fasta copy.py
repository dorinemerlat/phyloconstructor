#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def read_fasta(path):
    records = []

    header = None
    seq = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq)))

    return records


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def find_sequence(row, fasta_records, fasta_file):
    busco_id = str(row["busco_id"])
    sequence = str(row["sequence"])
    source = str(row["source"])

    busco_candidates = []
    exact_matches = []

    for header, seq in fasta_records:

        if busco_id not in header:
            continue

        busco_candidates.append(header)

        if sequence in header:
            exact_matches.append((header, seq))

    # Match unique => parfait
    if len(exact_matches) == 1:
        return exact_matches[0]

    # Plusieurs matches => problème
    if len(exact_matches) > 1:
        print(
            f"WARNING: multiple exact matches "
            f"for {row['specie_name']} "
            f"{source} "
            f"{row['data_id']} "
            f"{busco_id} "
            f"{sequence}"
        )

        print(f"FASTA: {fasta_file}")

        for h, _ in exact_matches:
            print(f"  MATCH: {h}")

        return exact_matches[0]

    # Cas genome
    if source == "ncbi_assemblies":

        if len(busco_candidates) == 1:
            print(
                f"WARNING: genome fallback "
                f"{row['specie_name']} "
                f"{row['data_id']} "
                f"{busco_id}"
            )

            return next(
                (h, s)
                for h, s in fasta_records
                if h == busco_candidates[0]
            )

        if len(busco_candidates) > 1:
            print(
                f"WARNING: multiple BUSCO candidates "
                f"for genome mode "
                f"{row['specie_name']} "
                f"{row['data_id']} "
                f"{busco_id}"
            )

            print(f"FASTA: {fasta_file}")

            for h in busco_candidates:
                print(f"  BUSCO: {h}")

    # Rien trouvé
    print(
        f"WARNING: sequence not found\n"
        f"  specie   = {row['specie_name']}\n"
        f"  source   = {source}\n"
        f"  data_id  = {row['data_id']}\n"
        f"  busco_id = {busco_id}\n"
        f"  sequence = {sequence}\n"
        f"  fasta    = {fasta_file}\n"
        f"  candidates_with_busco = {len(busco_candidates)}"
    )

    if len(busco_candidates) > 0:
        print("  Candidate headers:")

        for h in busco_candidates[:20]:
            print(f"    {h}")

        if len(busco_candidates) > 20:
            print(f"    ... {len(busco_candidates)-20} more")

    return None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--busco-sequences", nargs="+", required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.table, sep="\t")

    required = {
        "specie_taxid",
        "specie_name",
        "source",
        "data_id",
        "busco_id",
        "sequence",
    }

    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    fasta_files = []
    for item in args.busco_sequences:
        p = Path(item)
        if p.is_dir():
            fasta_files.extend(p.glob("**/*.fasta"))
            fasta_files.extend(p.glob("**/*.fa"))
            fasta_files.extend(p.glob("**/*.faa"))
        else:
            fasta_files.append(p)

    fasta_index = {}

    for fasta in fasta_files:
        name = fasta.name

        for _, row in table[["specie_name", "source", "data_id"]].drop_duplicates().iterrows():
            specie = str(row["specie_name"])
            source = str(row["source"])
            data_id = str(row["data_id"])

            if specie in name and source in name and data_id in name:
                key = (specie, source, data_id)
                fasta_index.setdefault(key, []).append(fasta)

    print(f"Rows in table: {len(table)}")
    print(f"Unique BUSCOs: {table['busco_id'].nunique()}")
    print(f"FASTA files loaded: {len(fasta_files)}")

    for busco_id, sub in table.groupby("busco_id"):
        output = f"{args.output_prefix}_{busco_id}.fasta"

        with open(output, "w") as out:
            for _, row in sub.iterrows():
                specie = str(row["specie_name"])
                taxid = str(row["specie_taxid"])
                source = str(row["source"])
                data_id = str(row["data_id"])
                sequence = str(row["sequence"])

                key = (specie, source, data_id)
                files = fasta_index.get(key, [])

                if not files:
                    print(f"WARNING: no BUSCO fasta for {key}")
                    continue

                found_header = None
                found_seq = None

                for fasta in files:
                    records = read_fasta(fasta)
                    found_header, found_seq = find_sequence(
                        row,
                        fasta_cache[fasta],
                        fasta
                    )

                    if found_seq is not None:
                        break

                if found_seq is None:
                    print(
                        f"WARNING: sequence not found: "
                        f"{specie} {source} {data_id} {busco_id} {sequence}"
                    )
                    continue

                clean_header = f"{specie}|{taxid}|{source}|{data_id}|{busco_id}|{sequence}"
                out.write(f">{clean_header}\n{wrap(found_seq)}\n")


if __name__ == "__main__":
    main()