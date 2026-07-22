process SELECT_BEST_BUSCO_HITS {
    tag "${label}/${specie}"
    memory '2 GB'
    time '1h'

    input:
    tuple val(label), val(taxid), val(specie), path(full_tables)

    output:
    tuple val(label), val(taxid), val(specie), path("${specie}_complete_table.tsv")

    script:
    """
    # Select the best BUSCO record for each ortholog across candidate datasets.
    select_best_busco_hits.py \\
        --full-tables ${full_tables.join(' ')} \\
        --output "${specie}_complete_table.tsv"
    """

    stub:
    """
    command -v select_best_busco_hits.py >/dev/null

    touch "${specie}_complete_table.tsv"
    """
}