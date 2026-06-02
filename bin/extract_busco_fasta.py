#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


HEADER_PROTEIN_MODE = re.compile(
    r"^(?P<specie>[^.]+)"
    r"\.(?P<source>[^.]+)"
    r"\.(?P<data_id>[^.]+)"
    r"\.(?P<orthogroup>[^.]+)"
    r"\.(?P<sequence_id>\S+)"
)

HEADER_ASSEMBLY_MODE = re.compile(
    r"^(?P<specie>.+?)"
    r"\.(?P<source>ncbi_assemblies)"
    r"\.(?P<data_id>GC[AF]_\d+\.\d+)"
    r"\.(?P<orthogroup>\d+at\d+)"
    r"_.*?\|(?P<sequence_id>[^:]+):"
)

HEADER_NCBI_PROTEOME_MODE = re.compile(
    r"^(?P<specie>.+?)"
    r"\.(?P<source>ncbi_proteomes)"
    r"\.(?P<data_id>GC[AF]_\d+\.\d+)"
    r"\.(?P<orthogroup>\d+at\d+)"
    r"\.(?P<sequence_id>\S+)"
)

def add_record(records, specie, source, data_id, orthogroup, sequence_id, seq):
    records.setdefault(specie, {})
    records[specie].setdefault(source, {})
    records[specie][source].setdefault(data_id, {})
    records[specie][source][data_id].setdefault(orthogroup, {})
    records[specie][source][data_id][orthogroup][sequence_id] = seq


def parse_fasta(fasta_file, records):
    for record in SeqIO.parse(fasta_file, "fasta"):
        match_assembly = HEADER_ASSEMBLY_MODE.match(record.description)
        match_ncbi_proteome = HEADER_NCBI_PROTEOME_MODE.match(record.description)
        match_protein = HEADER_PROTEIN_MODE.match(record.description)

        if match_assembly:
            match = match_assembly

        elif match_ncbi_proteome:
            match = match_ncbi_proteome

        elif match_protein:
            match = match_protein
            
        else:
            print(f"WARNING: header did not match expected format: {record.description}")
            exit(1)

        specie = match.group("specie")
        source = match.group("source")
        data_id = match.group("data_id")
        orthogroup = match.group("orthogroup")
        sequence_id = match.group("sequence_id")

        add_record(
            records=records,
            specie=specie,
            source=source,
            data_id=data_id,
            orthogroup=orthogroup,
            sequence_id=sequence_id,
            seq=record.seq,
        )
        
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--busco-sequences", nargs="+", required=True)
    parser.add_argument(
        "--output-dir",
        default="outgroups",
        help="Output directory for BUSCO FASTA files (default: outgroups)"
    )
    args = parser.parse_args()

    table = pd.read_csv(args.table, sep="\t")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Parsing BUSCO sequences...")
    single_copy_sequences = dict()
    multi_copy_sequences = dict()

    for file in args.busco_sequences:
        if "single_copy" in file:
            single_copy_sequences = parse_fasta(file, single_copy_sequences)

        elif "multi_copy" in file:
            multi_copy_sequences = parse_fasta(file, multi_copy_sequences)

    print("Extracting BUSCO sequences...")
    for _, row in table.iterrows():
        specie_taxid = str(row["specie_taxid"])
        specie_name = str(row["specie_name"])
        source = str(row["source"])
        data_id = str(row["data_id"])
        busco_id = str(row["busco_id"])
        status = str(row["status"])
        sequence_id = str(row["sequence"])

        if status == "Duplicated":
            fasta_record = (
                multi_copy_sequences
                .get(specie_name, {})
                .get(source, {})
                .get(data_id, {})
                .get(busco_id, {})
                .get(sequence_id)
            )
        else:
            fasta_record = (
                single_copy_sequences
                .get(specie_name, {})
                .get(source, {})
                .get(data_id, {})
                .get(busco_id, {})
                .get(sequence_id)
            )

        if fasta_record is None:
            print(row)
            print(
                f"WARNING: no FASTA records found\n"
                f"  specie      = {specie_name}\n"
                f"  source      = {source}\n"
                f"  data_id     = {data_id}\n"
                f"  busco_id    = {busco_id}\n"
                f"  sequence_id = {sequence_id}\n"
                f"  status      = {status}\n"
            )
            exit(1)

        output_file = output_dir / f"{args.output_prefix}_{busco_id}.fasta"

        new_record = SeqRecord(
            fasta_record,
            id=f"{specie_name} source={source} data_id={data_id} busco_id={busco_id} sequence_id={sequence_id}",
            description=""
        )

        with open(output_file, "a") as handle:
            SeqIO.write(new_record, handle, "fasta")


if __name__ == "__main__":
    main()