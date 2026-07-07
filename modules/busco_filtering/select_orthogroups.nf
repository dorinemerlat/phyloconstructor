process SELECT_ORTHOGROUPS {
    tag "$label"

    input:
    tuple val(label), path(table), val(dataset_size)

    output:
    tuple val(label), path("orthogroups_busco-c_*_gene-occupancy_*.tsv"), emit: tables
    tuple val(label), path("orthogroups_summary.tsv")
    tuple val(label), path("orthogroups_heatmap.{pdf,png,svg}")
    tuple val(label), path("orthogroups_heatmap_zoom_50_100.{pdf,png,svg}")
    tuple val(label), path("orthogroups_tradeoff.{pdf,png,svg}")

    script:
    """
    select_orthogroups.py \\
        --full-tables ${table.join(' ')} \\
        --output orthogroups \\
        --dataset-size ${dataset_size}
    """
}