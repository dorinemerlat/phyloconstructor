#!/usr/bin/env python3

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from ete3 import (
        Tree,
        TreeStyle,
        NodeStyle,
        TextFace,
        PieChartFace,
        faces,
    )
except ImportError:
    sys.exit(
        "\nERROR: ETE3 is required.\n\n"
        "Install it with:\n"
        "    pip install ete3 PyQt5\n"
    )


NON_TAXONOMY_COLUMNS = {
    "tree_label",
    "query_name",
    "specie",
    "species",
    "taxid",
}


QUARTET_COLORS = [
    "#377eb8",  # q1: topology shown
    "#e41a1c",  # q2: alternative topology 1
    "#4daf4a",  # q3: alternative topology 2
]

PIE_SIZE = 10


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate taxonomic clades and phylogenetic support on a Newick tree."
    )

    parser.add_argument("-i", "--input-tree", required=True, help="Input Newick tree.")
    parser.add_argument("-t", "--taxonomy", required=True, help="Taxonomy TSV table.")
    parser.add_argument("-o", "--output-prefix", default="taxonomic_tree", help="Output prefix. Default: taxonomic_tree")
    parser.add_argument("--tree-type", required=True, choices=["astral", "iqtree"], help="Type of input tree: astral or iqtree.")
    parser.add_argument("--outgroup-clade", help="Monophyletic taxonomic clade to use as outgroup.")
    parser.add_argument("--font-size", type=int, default=8, help="Leaf name font size. Default: 8")
    parser.add_argument("--clade-font-size", type=int, default=7, help="Clade and rank font size. Default: 7")
    parser.add_argument("--support-font-size", type=int, default=6, help="IQ-TREE bootstrap font size. Default: 6")
    parser.add_argument("--line-width", type=int, default=1, help="Branch line width. Default: 1")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution. Default: 300")

    return parser.parse_args()


def read_taxonomy_table(filename):
    """
    Read taxonomy TSV table.

    All columns except metadata columns are considered taxonomic ranks.
    """

    with open(filename, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header found in taxonomy table: {filename}")

        fieldnames = [field.strip() for field in reader.fieldnames]

        if "specie" in fieldnames:
            species_column = "specie"
        elif "species" in fieldnames:
            species_column = "species"
        else:
            raise ValueError(
                "The taxonomy table must contain a 'specie' or 'species' column."
            )

        taxonomy_columns = [
            column
            for column in fieldnames
            if column not in NON_TAXONOMY_COLUMNS
        ]

        rows = []

        for row in reader:
            cleaned_row = {
                key.strip(): (value or "").strip()
                for key, value in row.items()
                if key is not None
            }

            if cleaned_row.get(species_column):
                rows.append(cleaned_row)

    return rows, taxonomy_columns, species_column


def read_tree(filename):
    """
    Read a Newick tree with branch lengths and internal node names.
    """

    try:
        return Tree(filename, format=1, quoted_node_names=True)
    except Exception as error:
        raise RuntimeError(
            f"Could not parse Newick tree: {filename}\n\n"
            f"ETE3 error:\n{error}"
        )


def parse_astral_support(tree):
    """
    Extract q1, q2 and q3 from ASTRAL internal node annotations.

    Example:

        [q1=0.85;q2=0.10;q3=0.05;...]

    Values are stored as:
        node.q1
        node.q2
        node.q3
    """

    pattern = re.compile(
        r"q1=([0-9.eE+-]+);"
        r"q2=([0-9.eE+-]+);"
        r"q3=([0-9.eE+-]+)"
    )

    n_supported_nodes = 0

    for node in tree.traverse():
        if node.is_leaf() or not node.name:
            continue

        match = pattern.search(node.name)

        if not match:
            continue

        node.add_features(
            q1=float(match.group(1)),
            q2=float(match.group(2)),
            q3=float(match.group(3)),
        )

        n_supported_nodes += 1

    print(
        f"ASTRAL quartet support found for "
        f"{n_supported_nodes} internal nodes"
    )


def parse_iqtree_support(tree):
    """
    Extract bootstrap support from IQ-TREE internal node labels.

    Supported formats:

        98
            -> bootstrap = 98

        99.6/78
            -> bootstrap = 78

        100/100
            -> bootstrap = 100

    For labels containing two values separated by '/',
    the second value is used as bootstrap support.
    """

    pattern = re.compile(
        r"^\s*"
        r"([0-9]+(?:\.[0-9]+)?)"
        r"(?:/"
        r"([0-9]+(?:\.[0-9]+)?)"
        r")?"
        r"\s*$"
    )

    n_supported_nodes = 0

    for node in tree.traverse():
        if node.is_leaf() or not node.name:
            continue

        match = pattern.match(node.name)

        if not match:
            continue

        first_support = float(match.group(1))
        second_support = match.group(2)

        if second_support is not None:
            bootstrap = float(second_support)
        else:
            bootstrap = first_support

        node.add_feature(
            "bootstrap_support",
            bootstrap,
        )

        n_supported_nodes += 1

    print(
        f"IQ-TREE bootstrap support found for "
        f"{n_supported_nodes} internal nodes"
    )


def parse_phylogenetic_support(tree, tree_type):
    """
    Parse node support according to tree type.
    """

    if tree_type == "astral":
        parse_astral_support(tree)

    elif tree_type == "iqtree":
        parse_iqtree_support(tree)


def add_astral_quartet_pie(node):
    """
    Add a small q1/q2/q3 pie chart to an ASTRAL internal node.
    """

    if not all(
        hasattr(node, attribute)
        for attribute in ("q1", "q2", "q3")
    ):
        return

    values = [
        node.q1 * 100,
        node.q2 * 100,
        node.q3 * 100,
    ]

    pie_face = PieChartFace(
        values,
        width=PIE_SIZE,
        height=PIE_SIZE,
        colors=QUARTET_COLORS,
    )

    faces.add_face_to_node(
        pie_face,
        node,
        column=0,
        position="branch-top",
    )


def format_support_value(value):
    """
    Format support values without unnecessary '.0'.

    Examples:
        100.0 -> 100
        78.0  -> 78
        99.6  -> 99.6
    """

    if value.is_integer():
        return str(int(value))

    return f"{value:g}"


def add_iqtree_bootstrap(node, support_font_size):
    """
    Display IQ-TREE bootstrap support.

    Bootstrap == 100:
        show a black circle directly on the node.

    Bootstrap < 100:
        show the numerical bootstrap value above the branch.
    """

    if not hasattr(node, "bootstrap_support"):
        return

    if node.bootstrap_support == 100:
        node.img_style["size"] = 6
        node.img_style["shape"] = "circle"
        node.img_style["fgcolor"] = "black"

    else:
        support_face = TextFace(
            format_support_value(node.bootstrap_support),
            fsize=support_font_size,
            ftype="Arial",
            fgcolor="black",
        )

        support_face.margin_right = 2
        support_face.margin_bottom = 1

        faces.add_face_to_node(
            support_face,
            node,
            column=0,
            position="branch-top",
        )

def add_phylogenetic_support(
    node,
    tree_type,
    support_font_size,
):
    """
    Add graphical support according to tree type.
    """

    if node.is_leaf():
        return

    if tree_type == "astral":
        add_astral_quartet_pie(node)

    elif tree_type == "iqtree":
        add_iqtree_bootstrap(
            node,
            support_font_size,
        )


def build_taxonomic_groups(
    taxonomy_rows,
    taxonomy_columns,
    species_column,
    tree,
):
    """
    Build:

        groups[rank][clade] = set(species)
    """

    tree_species = set(tree.get_leaf_names())
    groups = defaultdict(lambda: defaultdict(set))

    for row in taxonomy_rows:
        species = row[species_column]

        if species not in tree_species:
            continue

        for rank in taxonomy_columns:
            clade = row.get(rank, "").strip()

            if clade:
                groups[rank][clade].add(species)

    return groups


def test_monophyly(tree, target_species):
    """
    Test strict monophyly.

    A group is monophyletic only if:

        leaves(MRCA(group)) == group
    """

    target_species = set(target_species)
    node = tree.get_common_ancestor(list(target_species))
    descendant_species = set(node.get_leaf_names())

    extra_species = descendant_species - target_species
    missing_species = target_species - descendant_species

    is_monophyletic = not extra_species and not missing_species

    return is_monophyletic, node, extra_species, missing_species


def root_tree_on_clade(tree, groups, root_clade):
    """
    Root the tree on the branch leading to a monophyletic clade.
    """

    matches = []

    for rank, rank_groups in groups.items():
        if root_clade in rank_groups:
            matches.append(
                (
                    rank,
                    rank_groups[root_clade],
                )
            )

    if not matches:
        raise ValueError(
            f"Outgroup clade '{root_clade}' "
            "was not found in the taxonomy table."
        )

    rank, species_set = max(
        matches,
        key=lambda x: len(x[1]),
    )

    if len(species_set) < 2:
        raise ValueError(
            f"Outgroup clade '{root_clade}' "
            "contains fewer than two species."
        )

    (
        is_monophyletic,
        node,
        extra_species,
        missing_species,
    ) = test_monophyly(
        tree,
        species_set,
    )

    if not is_monophyletic:
        raise ValueError(
            f"Outgroup clade '{root_clade}' ({rank}) "
            "is not monophyletic.\n"
            f"Extra species in MRCA: "
            f"{', '.join(sorted(extra_species))}"
        )

    print(
        f"Rooting tree on clade: {root_clade} "
        f"({rank}, {len(species_set)} species)"
    )

    tree.set_outgroup(node)


def find_taxonomic_clades(
    tree,
    groups,
    taxonomy_columns,
):
    """
    Test all taxonomic groups represented by at least two species.

    Monophyletic groups are assigned to their exact MRCA.

    Non-monophyletic groups are also assigned to their MRCA.
    """

    node_annotations = defaultdict(list)
    report_rows = []

    for rank in taxonomy_columns:

        for clade, species_set in groups.get(
            rank,
            {},
        ).items():

            if len(species_set) < 2:
                continue

            (
                is_monophyletic,
                node,
                extra_species,
                missing_species,
            ) = test_monophyly(
                tree,
                species_set,
            )

            node_annotations[node].append({
                "rank": rank,
                "clade": clade,
                "monophyletic": is_monophyletic,
            })

            report_rows.append({
                "rank": rank,
                "clade": clade,
                "n_species": len(species_set),
                "monophyletic": (
                    "yes"
                    if is_monophyletic
                    else "no"
                ),
                "extra_species": ";".join(
                    sorted(extra_species)
                ),
                "missing_species": ";".join(
                    sorted(missing_species)
                ),
            })

            if not is_monophyletic:
                print(
                    f"WARNING: {clade} ({rank}) "
                    "is not monophyletic"
                )

                if extra_species:
                    print(
                        "    Extra species in MRCA: "
                        + ", ".join(
                            sorted(extra_species)
                        )
                    )

    return node_annotations, report_rows


def merge_annotations(
    node_annotations,
    taxonomy_columns,
):
    """
    Merge taxonomic levels assigned to the same node.

    Monophyletic and non-monophyletic annotations are kept separate.
    """

    rank_order = {
        rank: index
        for index, rank in enumerate(
            taxonomy_columns
        )
    }

    merged = {}

    for node, annotations in node_annotations.items():

        annotations = sorted(
            annotations,
            key=lambda item: rank_order.get(
                item["rank"],
                10**9,
            ),
        )

        mono_clades = []
        mono_ranks = []
        nonmono_clades = []
        nonmono_ranks = []

        for annotation in annotations:

            clade = annotation["clade"]
            rank = annotation["rank"]

            if annotation["monophyletic"]:

                if clade not in mono_clades:
                    mono_clades.append(clade)

                if rank not in mono_ranks:
                    mono_ranks.append(rank)

            else:

                if clade not in nonmono_clades:
                    nonmono_clades.append(clade)

                if rank not in nonmono_ranks:
                    nonmono_ranks.append(rank)

        merged[node] = {
            "mono_clade_label": " / ".join(
                mono_clades
            ),
            "mono_rank_label": " / ".join(
                mono_ranks
            ),
            "nonmono_clade_label": " / ".join(
                nonmono_clades
            ),
            "nonmono_rank_label": " / ".join(
                nonmono_ranks
            ),
        }

    return merged


def write_monophyly_report(
    filename,
    report_rows,
):
    """
    Write monophyly report.
    """

    fieldnames = [
        "rank",
        "clade",
        "n_species",
        "monophyletic",
        "extra_species",
        "missing_species",
    ]

    with open(
        filename,
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(report_rows)


def add_annotations_to_newick_nodes(
    merged_annotations,
):
    """
    Add taxonomic annotations to internal node names.
    """

    for node, annotation in merged_annotations.items():

        annotations = []

        if annotation["mono_clade_label"]:

            annotations.append(
                "MONOPHYLETIC_CLADE="
                f"{annotation['mono_clade_label']}"
            )

            annotations.append(
                "MONOPHYLETIC_RANK="
                f"{annotation['mono_rank_label']}"
            )

        if annotation["nonmono_clade_label"]:

            annotations.append(
                "NON_MONOPHYLETIC_CLADE="
                f"{annotation['nonmono_clade_label']}"
            )

            annotations.append(
                "NON_MONOPHYLETIC_RANK="
                f"{annotation['nonmono_rank_label']}"
            )

        taxonomy_annotation = " | ".join(
            annotations
        )

        if not taxonomy_annotation:
            continue

        if node.name:
            node.name = (
                f"{node.name.strip()} | "
                f"{taxonomy_annotation}"
            )
        else:
            node.name = taxonomy_annotation


def write_annotated_tree(tree, filename):
    """
    Write annotated Newick tree.
    """

    tree.write(
        outfile=filename,
        format=1,
        quoted_node_names=True,
        format_root_node=True,
    )


def set_branch_style(
    tree,
    line_width,
):
    """
    Apply common branch style.
    """

    for node in tree.traverse():

        node_style = NodeStyle()

        node_style["size"] = 0
        node_style["hz_line_width"] = line_width
        node_style["vt_line_width"] = line_width

        node.set_style(
            node_style
        )


def style_tree(
    tree,
    merged_annotations,
    tree_type,
    font_size,
    clade_font_size,
    support_font_size,
    line_width,
):
    """
    Tree with taxonomic annotations and phylogenetic support.
    """

    tree_style = TreeStyle()

    tree_style.mode = "r"
    tree_style.show_leaf_name = False
    tree_style.show_scale = True
    tree_style.branch_vertical_margin = 4

    set_branch_style(
        tree,
        line_width,
    )

    def layout(node):

        # ASTRAL quartet pie or IQ-TREE bootstrap.
        add_phylogenetic_support(
            node,
            tree_type,
            support_font_size,
        )

        # Leaf names.
        if node.is_leaf():

            leaf_face = TextFace(
                node.name,
                fsize=font_size,
                ftype="Arial",
                fgcolor="black",
                fstyle="italic",
            )

            leaf_face.margin_left = 3

            faces.add_face_to_node(
                leaf_face,
                node,
                column=0,
                position="branch-right",
            )

        if node not in merged_annotations:
            return

        annotation = merged_annotations[node]

        # Monophyletic clades.
        if annotation["mono_clade_label"]:

            clade_face = TextFace(
                annotation[
                    "mono_clade_label"
                ],
                fsize=clade_font_size,
                ftype="Arial",
                fgcolor="#377eb8",
            )

            clade_face.margin_left = 2
            clade_face.margin_bottom = 1

            faces.add_face_to_node(
                clade_face,
                node,
                column=1,
                position="branch-top",
            )

            rank_face = TextFace(
                annotation[
                    "mono_rank_label"
                ],
                fsize=clade_font_size,
                ftype="Arial",
                fgcolor="#e34a33",
            )

            rank_face.margin_left = 2
            rank_face.margin_top = 1

            faces.add_face_to_node(
                rank_face,
                node,
                column=1,
                position="branch-bottom",
            )

        # Non-monophyletic clades.
        if annotation[
            "nonmono_clade_label"
        ]:

            nonmono_face = TextFace(
                annotation[
                    "nonmono_clade_label"
                ],
                fsize=clade_font_size,
                ftype="Arial",
                fgcolor="#e67e22",
            )

            nonmono_face.margin_left = 4
            nonmono_face.margin_bottom = 1

            faces.add_face_to_node(
                nonmono_face,
                node,
                column=2,
                position="branch-top",
            )

            nonmono_rank_face = TextFace(
                annotation[
                    "nonmono_rank_label"
                ],
                fsize=clade_font_size,
                ftype="Arial",
                fgcolor="#a04000",
            )

            nonmono_rank_face.margin_left = 4
            nonmono_rank_face.margin_top = 1

            faces.add_face_to_node(
                nonmono_rank_face,
                node,
                column=2,
                position="branch-bottom",
            )

    tree_style.layout_fn = layout

    return tree_style

def get_tree_scale(tree, target_width=1200):
    """
    Calculate a scale so that the longest root-to-tip distance
    occupies approximately target_width pixels.
    """

    max_distance = max(
        tree.get_distance(leaf)
        for leaf in tree.iter_leaves()
    )

    if max_distance <= 0:
        return 100

    return target_width / max_distance

def style_tree_without_annotations(
    tree,
    tree_type,
    font_size,
    support_font_size,
    line_width,
):
    """
    Tree without taxonomic clade annotations,
    but with phylogenetic support.
    """

    tree_style = TreeStyle()

    tree_style.mode = "r"
    tree_style.show_leaf_name = False
    tree_style.show_scale = True
    tree_style.branch_vertical_margin = 4

    tree_style.scale = get_tree_scale(tree, target_width=1000)

    tree_style.margin_top = 20
    tree_style.margin_bottom = 40
    tree_style.margin_left = 20
    tree_style.margin_right = 20

    set_branch_style(
        tree,
        line_width,
    )

    def layout(node):

        # ASTRAL quartet pie or IQ-TREE bootstrap.
        add_phylogenetic_support(
            node,
            tree_type,
            support_font_size,
        )

        if node.is_leaf():

            leaf_face = TextFace(
                node.name,
                fsize=font_size,
                ftype="Arial",
                fgcolor="black",
                fstyle="italic",
            )

            leaf_face.margin_left = 3

            faces.add_face_to_node(
                leaf_face,
                node,
                column=0,
                position="branch-right",
            )

    tree_style.layout_fn = layout

    return tree_style


def render_tree(
    tree,
    tree_style,
    output_prefix,
    width_mm,
    dpi,
):
    """
    Render PNG, SVG and PDF.

    If width_mm is None, the output size is determined automatically
    from the tree geometry and TreeStyle.scale.
    """

    outputs = {
        "PNG": f"{output_prefix}.png",
        "SVG": f"{output_prefix}.svg",
        "PDF": f"{output_prefix}.pdf",
    }

    for filetype, filename in outputs.items():
        print(f"Rendering {filetype}: {filename}")

        kwargs = {
            "tree_style": tree_style,
        }

        # Only force an output width when explicitly requested.
        if width_mm is not None:
            kwargs["w"] = width_mm
            kwargs["units"] = "mm"

        if filetype == "PNG":
            kwargs["dpi"] = dpi

        tree.render(
            filename,
            **kwargs,
        )


def main():
    args = parse_args()

    print(
        f"Reading {args.tree_type.upper()} tree: "
        f"{args.input_tree}"
    )

    tree = read_tree(
        args.input_tree
    )

    print(
        f"Tree contains {len(tree)} leaves"
    )

    # Parse ASTRAL or IQ-TREE support before
    # modifying internal node names.
    parse_phylogenetic_support(
        tree,
        args.tree_type,
    )

    print(
        f"\nReading taxonomy table: "
        f"{args.taxonomy}"
    )

    (
        taxonomy_rows,
        taxonomy_columns,
        species_column,
    ) = read_taxonomy_table(
        args.taxonomy
    )

    print(
        f"Taxonomy table contains "
        f"{len(taxonomy_rows)} species"
    )

    print(
        f"Taxonomic levels detected: "
        f"{len(taxonomy_columns)}"
    )

    groups = build_taxonomic_groups(
        taxonomy_rows,
        taxonomy_columns,
        species_column,
        tree,
    )

    if args.outgroup_clade:

        root_tree_on_clade(
            tree,
            groups,
            args.outgroup_clade,
        )

    total_groups = sum(
        len(groups[rank])
        for rank in taxonomy_columns
    )

    print(
        f"\nTaxonomic groups detected: "
        f"{total_groups}"
    )

    print(
        "Testing monophyly..."
    )

    (
        node_annotations,
        report_rows,
    ) = find_taxonomic_clades(
        tree,
        groups,
        taxonomy_columns,
    )

    merged_annotations = merge_annotations(
        node_annotations,
        taxonomy_columns,
    )

    monophyletic_groups = sum(
        row["monophyletic"] == "yes"
        for row in report_rows
    )

    non_monophyletic_groups = sum(
        row["monophyletic"] == "no"
        for row in report_rows
    )

    print(
        "\nMonophyly summary:"
    )

    print(
        f"    Monophyletic groups     : "
        f"{monophyletic_groups}"
    )

    print(
        f"    Non-monophyletic groups : "
        f"{non_monophyletic_groups}"
    )

    print(
        f"    Annotated nodes         : "
        f"{len(merged_annotations)}"
    )

    report_file = (
        f"{args.output_prefix}"
        "_monophyly_report.tsv"
    )

    write_monophyly_report(
        report_file,
        report_rows,
    )

    print(
        f"\nWritten report: "
        f"{report_file}"
    )

    # --------------------------------------------------------
    # Render tree with clade annotations
    # --------------------------------------------------------

    tree_style = style_tree(
        tree,
        merged_annotations,
        tree_type=args.tree_type,
        font_size=args.font_size,
        clade_font_size=args.clade_font_size,
        support_font_size=args.support_font_size,
        line_width=args.line_width,
    )

    render_tree(
        tree,
        tree_style,
        args.output_prefix,
        300,
        args.dpi,
    )

    # --------------------------------------------------------
    # Save and render tree without clade annotations
    # --------------------------------------------------------

    plain_prefix = (
        f"{args.output_prefix}"
        "_without_clade_annotations"
    )

    plain_tree_file = (
        f"{plain_prefix}.nwk"
    )

    tree.write(
        outfile=plain_tree_file,
        format=1,
        quoted_node_names=True,
        format_root_node=True,
    )

    print(
        "Written tree without clade annotations: "
        f"{plain_tree_file}"
    )

    plain_tree_style = style_tree_without_annotations(
        tree,
        tree_type=args.tree_type,
        font_size=args.font_size,
        support_font_size=args.support_font_size,
        line_width=args.line_width,
    )

    render_tree(
        tree,
        plain_tree_style,
        plain_prefix,
        500,
        args.dpi,
    )

    # --------------------------------------------------------
    # Add taxonomic annotations to Newick nodes
    # --------------------------------------------------------

    add_annotations_to_newick_nodes(
        merged_annotations
    )

    output_tree = (
        f"{args.output_prefix}"
        "_annotated.nwk"
    )

    write_annotated_tree(
        tree,
        output_tree,
    )

    print(
        f"Written annotated tree: "
        f"{output_tree}"
    )

    print(
        "Done."
    )


if __name__ == "__main__":
    main()