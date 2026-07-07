process SELECT_BEST_BUSCO_HITS {
    tag "${label}/${specie}"

    input:
    tuple val(label), val(taxid), val(specie), path(full_tables)

    output:
    tuple val(label), val(taxid), val(specie), path("${specie}_complete_table.tsv")

    script:
    """
    select_best_busco_hits.py \\
        --full-tables ${full_tables.join(' ')} \\
        --output ${specie}_complete_table.tsv
    """
}