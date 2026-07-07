#!/usr/bin/env python3

import argparse
import csv
import re


EXPECTED_BUSCO_COLUMNS = [
    "busco_id",
    "status",
    "sequence",
    "score",
    "length",
    "orthodb_url",
    "description",
]


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def clean_colname(name):
    return re.sub(r"_+", "_", name.strip().lower().replace(" ", "_"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-table", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--specie", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--data-id", required=True)
    args = parser.parse_args()

    header = None
    clean_header = None
    rows = []
    best_duplicated = {}

    with open(args.full_table) as handle:
        for line in handle:
            line = line.rstrip("\n")

            if line.startswith("# Busco id"):
                header = line.lstrip("# ").split("\t")
                clean_header = [clean_colname(col) for col in header]
                continue

            if line.startswith("#") or not line.strip():
                continue

            if header is None:
                raise ValueError(
                    "Could not find BUSCO header line starting with '# Busco id'"
                )

            cols = line.split("\t")
            row = dict(zip(clean_header, cols))

            status = row.get("status", "")
            if status not in {"Complete", "Duplicated"}:
                continue

            out_row = {
                "specie_taxid": args.taxid,
                "specie_name": args.specie,
                "source": args.source,
                "data_id": args.data_id,
            }

            for col in EXPECTED_BUSCO_COLUMNS:
                out_row[col] = row.get(col, "")

            if status == "Complete":
                rows.append(out_row)

            elif status == "Duplicated":
                busco_id = out_row.get("busco_id", "")
                current_score = to_float(out_row.get("score", ""))
                current_length = to_float(out_row.get("length", ""))

                previous = best_duplicated.get(busco_id)

                if previous is None:
                    best_duplicated[busco_id] = out_row
                    continue

                previous_score = to_float(previous.get("score", ""))
                previous_length = to_float(previous.get("length", ""))

                if current_score > previous_score:
                    best_duplicated[busco_id] = out_row
                elif current_score == previous_score and current_length > previous_length:
                    best_duplicated[busco_id] = out_row

    rows.extend(best_duplicated.values())

    out_fields = [
        "specie_taxid",
        "specie_name",
        "source",
        "data_id",
    ] + EXPECTED_BUSCO_COLUMNS

    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()