process SELECT_ORTHOGROUPS {
    tag ""

    input:
    tuple path(table), val(dataset_size)

    output:
    path("orthogroups_busco-c_*_gene-occupancy_*.tsv"), emit: tables
    path("orthogroups_summary.tsv")
    path("orthogroups_heatmap.{pdf,png,svg}")
    path("orthogroups_tradeoff.{pdf,png,svg}")

    script:
    """
    select_orthogroups.py \\
        --full-tables ${table.join(' ')} \\
        --output orthogroups \\
        --dataset-size ${dataset_size}
    """
}