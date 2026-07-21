#!/usr/bin/env python3

import argparse
import csv
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


# Order used for standard taxonomic ranks.
# Additional intermediate ranks such as:
#   phylum-1
#   phylum-2
#   class-1
# are inserted dynamically.
RANK_ORDER = [
    "superkingdom",
    "kingdom",
    "subkingdom",
    "superphylum",
    "phylum",
    "subphylum",
    "infraphylum",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "species",
    "subspecies",
    "varietas",
    "forma",
]


def fetch_url(url, timeout=30, max_retries=3):
    """
    Download a URL with retry handling.
    """

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "newick-taxonomy-fetcher/2.0"
        },
    )

    for attempt in range(1, max_retries + 1):

        try:
            with urllib.request.urlopen(
                req,
                timeout=timeout,
            ) as response:

                return response.read()

        except urllib.error.HTTPError as e:

            if (
                e.code in {429, 500, 502, 503, 504}
                and attempt < max_retries
            ):

                wait_time = attempt * 5

                print(
                    f"WARNING: HTTP {e.code}, "
                    f"retrying in {wait_time} seconds..."
                )

                time.sleep(wait_time)

            else:
                raise

        except urllib.error.URLError:

            if attempt < max_retries:

                wait_time = attempt * 5

                print(
                    "WARNING: network error, "
                    f"retrying in {wait_time} seconds..."
                )

                time.sleep(wait_time)

            else:
                raise


def convert_tree_label(label):
    """
    Convert a tree label into a taxonomic query name.

    Rules:
        --  -> ". "
        -   -> " "

    The first character is capitalized.

    Examples:
        scoterpes-sp--aums15047
            -> Scoterpes sp. aums15047

        julidae-sp--jj-2019
            -> Julidae sp. jj 2019

        drosophila-melanogaster
            -> Drosophila melanogaster
    """

    placeholder = "\0DOUBLE_HYPHEN\0"

    name = label.replace("--", placeholder)
    name = name.replace("-", " ")
    name = name.replace(placeholder, ". ")

    # Remove duplicated whitespace.
    name = re.sub(r"\s+", " ", name).strip()

    # Capitalize only the first character.
    if name:
        name = name[0].upper() + name[1:]

    return name


def extract_leaf_names_from_newick(newick):
    """
    Extract terminal leaf labels from a Newick string.

    A leaf label occurs immediately after:
        (
        ,

    This avoids extracting internal node labels and annotations.

    Supports:
        species-name
        'Species name'
    """

    pattern = re.compile(
        r"(?:\(|,)\s*"
        r"("
        r"'(?:[^']|'')*'"
        r"|"
        r"[^(),:;\s]+"
        r")"
    )

    labels = []

    for match in pattern.finditer(newick):

        label = match.group(1).strip()

        # Remove Newick quotes.
        if (
            label.startswith("'")
            and label.endswith("'")
        ):
            label = label[1:-1]
            label = label.replace("''", "'")

        if label:
            labels.append(label)

    # Remove duplicates while preserving order.
    return list(dict.fromkeys(labels))


def read_newick(newick_file):
    """
    Read the complete Newick tree as a string.
    """

    with open(
        newick_file,
        "r",
        encoding="utf-8",
    ) as fh:

        return fh.read().strip()


def resolve_name_to_taxid(name):
    """
    Resolve a scientific name to an NCBI Taxonomy taxid.

    First:
        exact Scientific Name search

    Then:
        broader taxonomy search
    """

    base_url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/esearch.fcgi"
    )

    queries = [
        f'"{name}"[Scientific Name]',
        name,
    ]

    for query in queries:

        params = urllib.parse.urlencode({
            "db": "taxonomy",
            "term": query,
            "retmode": "xml",
        })

        url = f"{base_url}?{params}"

        data = fetch_url(url)

        root = ET.fromstring(data)

        ids = [
            element.text
            for element in root.findall(".//IdList/Id")
            if element.text
        ]

        if ids:
            return ids[0]

    raise ValueError(
        f"No taxid found for name: {name}"
    )


def make_intermediate_rank(
    previous_rank,
    intermediate_counters,
):
    """
    Generate a synthetic intermediate rank name.

    Examples:
        previous official rank = phylum

        first no-rank taxon:
            phylum-1

        second no-rank taxon:
            phylum-2

    The counter is specific to the previous official rank.
    """

    intermediate_counters[previous_rank] = (
        intermediate_counters.get(previous_rank, 0)
        + 1
    )

    return (
        f"{previous_rank}-"
        f"{intermediate_counters[previous_rank]}"
    )


def parse_taxonomic_lineage(taxon):
    """
    Parse the complete NCBI lineage.

    Standard ranks are kept as provided by NCBI.

    Taxa with rank:
        no rank
        clade

    are converted to synthetic intermediate ranks based on
    the previous official rank.

    Example:

        phylum      Arthropoda
        no rank     Mandibulata
        clade       Myriapoda
        class       Diplopoda

    becomes:

        phylum      Arthropoda
        phylum-1    Mandibulata
        phylum-2    Myriapoda
        class       Diplopoda

    Returns:
        lineage: dict
        ordered_ranks: list
    """

    lineage = {}
    ordered_ranks = []

    previous_official_rank = None
    intermediate_counters = {}

    lineage_taxa = list(
        taxon.findall("./LineageEx/Taxon")
    )

    # Add current taxon at the end.
    lineage_taxa.append(taxon)

    for lineage_taxon in lineage_taxa:

        rank = lineage_taxon.findtext(
            "Rank",
            default="",
        ).strip().lower()

        name = lineage_taxon.findtext(
            "ScientificName",
            default="",
        ).strip()

        if not name:
            continue

        # Standard ranked taxon.
        if rank and rank not in {
            "no rank",
            "clade",
        }:

            output_rank = rank

            previous_official_rank = rank

            # Reset counter for this rank.
            intermediate_counters[rank] = 0

        else:

            # There is no previous official rank.
            # This can happen near the root of taxonomy.
            if previous_official_rank is None:
                previous_official_rank = "root"

            output_rank = make_intermediate_rank(
                previous_official_rank,
                intermediate_counters,
            )

        lineage[output_rank] = name
        ordered_ranks.append(output_rank)

    return lineage, ordered_ranks


def fetch_taxonomy(taxid):
    """
    Fetch complete taxonomy information from NCBI.
    """

    params = urllib.parse.urlencode({
        "db": "taxonomy",
        "id": taxid,
        "retmode": "xml",
    })

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi"
        f"?{params}"
    )

    data = fetch_url(url)

    root = ET.fromstring(data)

    taxon = root.find(".//Taxon")

    if taxon is None:
        raise ValueError(
            f"No taxonomy information found "
            f"for taxid: {taxid}"
        )

    scientific_name = taxon.findtext(
        "ScientificName",
        default="",
    ).strip()

    resolved_taxid = taxon.findtext(
        "TaxId",
        default=str(taxid),
    ).strip()

    lineage, ordered_ranks = (
        parse_taxonomic_lineage(taxon)
    )

    return {
        "specie": scientific_name,
        "taxid": resolved_taxid,
        "lineage": lineage,
        "ordered_ranks": ordered_ranks,
    }


def newick_quote_name(name):
    """
    Quote a taxon name for safe use in a Newick tree.

    Spaces in Newick labels can be problematic, so names are
    written between single quotes.

    Example:
        Drosophila melanogaster

    becomes:
        'Drosophila melanogaster'
    """

    # Escape single quotes according to Newick convention.
    escaped = name.replace("'", "''")

    return f"'{escaped}'"


def replace_leaf_names(newick, name_mapping):
    """
    Replace terminal leaf labels in a Newick tree.

    Only leaf labels immediately following:
        (
        ,

    are replaced.

    Internal node annotations and support values are untouched.

    name_mapping:
        original_label -> new_scientific_name
    """

    pattern = re.compile(
        r"(?P<prefix>\(|,)"
        r"(?P<space>\s*)"
        r"(?P<label>"
        r"'(?:[^']|'')*'"
        r"|"
        r"[^(),:;\s]+"
        r")"
    )

    def replacement(match):

        prefix = match.group("prefix")
        whitespace = match.group("space")
        raw_label = match.group("label")

        label = raw_label

        if (
            label.startswith("'")
            and label.endswith("'")
        ):

            label = label[1:-1]
            label = label.replace("''", "'")

        if label not in name_mapping:
            return match.group(0)

        new_name = name_mapping[label]

        return (
            prefix
            + whitespace
            + newick_quote_name(new_name)
        )

    return pattern.sub(
        replacement,
        newick,
    )


def split_rank_name(rank):
    """
    Split a rank into:
        base rank
        intermediate index

    Examples:
        phylum
            -> ("phylum", 0)

        phylum-1
            -> ("phylum", 1)

        phylum-2
            -> ("phylum", 2)
    """

    match = re.match(
        r"^(.*)-([0-9]+)$",
        rank,
    )

    if match:

        return (
            match.group(1),
            int(match.group(2)),
        )

    return rank, 0


def get_rank_sort_key(rank):
    """
    Return a sorting key for taxonomy table columns.

    Intermediate ranks are placed immediately after their
    corresponding standard rank.

    Example:
        phylum
        phylum-1
        phylum-2
        subphylum
        class
        class-1
        order
    """

    base_rank, intermediate_index = (
        split_rank_name(rank)
    )

    try:
        base_index = RANK_ORDER.index(base_rank)

    except ValueError:
        # Unknown standard ranks are placed after known ranks.
        base_index = len(RANK_ORDER)

    return (
        base_index,
        intermediate_index,
        rank,
    )


def determine_output_ranks(results):
    """
    Determine every taxonomy column required across all species.

    Includes:
        standard ranks
        dynamically generated intermediate ranks
    """

    all_ranks = set()

    for result in results:

        all_ranks.update(
            result["lineage"].keys()
        )

    return sorted(
        all_ranks,
        key=get_rank_sort_key,
    )


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract leaves from a Newick tree, retrieve their "
            "complete NCBI taxonomy, generate a taxonomy table "
            "including intermediate no-rank/clade levels, and "
            "write a renamed Newick tree."
        )
    )

    parser.add_argument("-i", "--input", required=True, help="Input Newick tree.")

    parser.add_argument("-o", "--output", default="tree_taxonomy.tsv",
        help=("Output taxonomy TSV file. Default: tree_taxonomy.tsv"))

    parser.add_argument("-t", "--output-tree", default="tree_renamed.nwk",
     help=("Output Newick tree with renamed leaves. Default: tree_renamed.nwk"))

    args = parser.parse_args()

    sleep_time = 0.3

    # --------------------------------------------------------
    # Read tree
    # --------------------------------------------------------

    print(f"Reading tree: {args.input}")

    newick = read_newick(args.input)

    tree_labels = (extract_leaf_names_from_newick( newick ))

    print(f"Found {len(tree_labels)} unique terminal taxa")

    # --------------------------------------------------------
    # Fetch taxonomy
    # --------------------------------------------------------

    results = []

    # Mapping used to rename tree leaves.
    name_mapping = {}

    MAX_SPECIES_RETRIES = 5

    for index, tree_label in enumerate(tree_labels, start=1):

        query_name = convert_tree_label(tree_label)

        success = False

        for attempt in range(MAX_SPECIES_RETRIES):

            try:

                # Resolve name.
                taxid = resolve_name_to_taxid(query_name)

                time.sleep(sleep_time)

                # Retrieve taxonomy.
                taxonomy = fetch_taxonomy(taxid)

                scientific_name = taxonomy["specie"]

                results.append({
                    "tree_label": tree_label,
                    "query_name": query_name,
                    "specie": scientific_name,
                    "taxid": taxonomy["taxid"],
                    "lineage": taxonomy["lineage"],
                    "ordered_ranks": taxonomy["ordered_ranks"],
                })

                # Rename tree leaf with official NCBI name.
                name_mapping[tree_label] = scientific_name

                print(
                    f"{index} / {len(tree_labels)} : "
                    f"{scientific_name} found"
                )

                success = True
                break

            except Exception:

                # Retry silently with exponential backoff.
                if attempt < MAX_SPECIES_RETRIES - 1:

                    wait_time = min(
                        120,
                        5 * (2 ** attempt),
                    )

                    time.sleep(wait_time)

        # Failed after all retries.
        if not success:

            fallback_name = query_name

            results.append({
                "tree_label": tree_label,
                "query_name": query_name,
                "specie": fallback_name,
                "taxid": "NA",
                "lineage": {},
                "ordered_ranks": [],
            })

            name_mapping[tree_label] = fallback_name

            print(
                f"{index} / {len(tree_labels)} : "
                f"{query_name} not found"
            )

        time.sleep(sleep_time)

    # --------------------------------------------------------
    # Determine all taxonomy columns
    # --------------------------------------------------------

    output_ranks = determine_output_ranks(results)

    # --------------------------------------------------------
    # Write taxonomy table
    # --------------------------------------------------------

    fieldnames = [
        "tree_label",
        "query_name",
        "specie",
        "taxid",
    ] + output_ranks

    with open(args.output, "w") as out:

        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")

        writer.writeheader()

        for result in results:

            row = {
                "tree_label": (result["tree_label"]),
                "query_name": (result["query_name"]),
                "specie": (result["specie"]),
                "taxid": (result["taxid"])
            }

            for rank in output_ranks:

                row[rank] = (result["lineage"].get(rank, ""))

            writer.writerow(row)

    print(f"Taxonomy table written to: {args.output}")

    # --------------------------------------------------------
    # Rename tree
    # --------------------------------------------------------

    renamed_newick = replace_leaf_names(newick, name_mapping)

    with open(args.output_tree, "w") as out:

        out.write(renamed_newick)

        if not renamed_newick.endswith("\n"):
            out.write("\n")

    print(f"Renamed tree written to: {args.output_tree}")

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    successful = sum(
        result["taxid"] != "NA"
        for result in results
    )

    failed = (len(results) - successful )

    print()
    print("========================================")
    print(f"Taxa in tree       : {len(results)}")
    print(f"Resolved by NCBI   : {successful}")
    print(f"Unresolved         : {failed}")
    print(f"Taxonomy table     : {args.output}")
    print(f"Renamed tree       : {args.output_tree}")
    print("========================================")


if __name__ == "__main__":
    main()