#!/usr/bin/env python3

import argparse
import sys

from ete3 import Tree


def parse_args():
    parser = argparse.ArgumentParser(description="Prune a Newick tree to keep only the species listed in a text file.")
    parser.add_argument("-i", "--input-tree", required=True, help="Input Newick tree.")
    parser.add_argument("-l", "--list", required=True, help="Text file containing one species name per line.")
    parser.add_argument("-o", "--output-tree", required=True, help="Output pruned Newick tree.")
    return parser.parse_args()


def read_species_list(filename):
    species = []

    with open(filename, "r", encoding="utf-8") as handle:
        for line in handle:
            name = line.strip()

            if not name or name.startswith("#"):
                continue

            species.append(name)

    return species


def main():
    args = parse_args()

    # format=1 preserves internal node names such as IQ-TREE support:
    # 100/100, 99.6/78, etc.
    tree = Tree(args.input_tree, format=1, quoted_node_names=True)

    species_to_keep = read_species_list(args.list)

    tree_species = set(tree.get_leaf_names())
    requested_species = set(species_to_keep)

    found_species = requested_species & tree_species
    missing_species = requested_species - tree_species

    print(f"Species in input tree : {len(tree_species)}")
    print(f"Species requested     : {len(requested_species)}")
    print(f"Species found         : {len(found_species)}")
    print(f"Species not found     : {len(missing_species)}")

    if missing_species:
        print("\nWARNING: species not found in the tree:")

        for species in sorted(missing_species):
            print(f"    {species}")

    if len(found_species) < 2:
        sys.exit("\nERROR: fewer than two requested species were found in the tree.")

    # Keep only the requested leaves.
    #
    # preserve_branch_length=True ensures that when intermediate nodes
    # disappear, their branch lengths are added together so that the
    # original distances between retained taxa are preserved.
    tree.prune(
        list(found_species),
        preserve_branch_length=True,
    )

    tree.write(
        outfile=args.output_tree,
        format=1,
        quoted_node_names=True,
        format_root_node=True,
    )

    print(f"\nSpecies retained: {len(tree)}")
    print(f"Written: {args.output_tree}")


if __name__ == "__main__":
    main()