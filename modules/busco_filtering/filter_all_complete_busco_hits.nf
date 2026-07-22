process FILTER_ALL_COMPLETE_BUSCO_HITS {
    tag "${specie}/${source}/${data_id}"
    memory '1 GB'
    time '1h'

    input:
    tuple val(taxid), val(specie), val(data_id), path(full_table), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}_complete_table.tsv"), val(source)

    script:
    """
    # Retain all complete BUSCO records, including duplicated BUSCOs.
    filter_all_complete_busco_hits.py \\
        --full-table "${full_table}" \\
        --output "${specie}_${source}_${data_id}_complete_table.tsv" \\
        --taxid "${taxid}" \\
        --specie "${specie}" \\
        --source "${source}" \\
        --data-id "${data_id}"
    """

    stub:
    """
    command -v filter_all_complete_busco_hits.py >/dev/null

    touch "${specie}_${source}_${data_id}_complete_table.tsv"
    """
}