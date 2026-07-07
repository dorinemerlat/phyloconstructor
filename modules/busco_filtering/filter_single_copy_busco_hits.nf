process FILTER_SINGLE_COPY_BUSCO_HITS {
    tag "${specie}"

    cpus 1
    memory '1 GB'

    input:
    tuple val(taxid), val(specie), val(data_id), path(full_table), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}_single_copy_table.tsv"), val(source)

    script:
    """
    filter_single_copy_busco_hits.py \\
        --full-table ${full_table} \\
        --output ${specie}_${source}_${data_id}_single_copy_table.tsv \\
        --taxid ${taxid} \\
        --specie ${specie} \\
        --source ${source} \\
        --data-id ${data_id}
    """
}