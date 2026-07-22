process SELECT_ORTHOGROUPS {
    tag "${label}"
    cpus 2
    memory '4 GB'
    time '2h'

    input:
    tuple val(label), path(table), val(dataset_size)

    output:
    tuple val(label), path("orthogroups_busco-c_*_gene-occupancy_*.tsv"), emit: tables
    tuple val(label), path("orthogroups_summary.tsv"), emit: summary
    tuple val(label), path("orthogroups_heatmap.{pdf,png,svg}"), emit: heatmap
    tuple val(label), path("orthogroups_heatmap_zoom_50_100.{pdf,png,svg}"), emit: heatmap_zoom_50_100
    tuple val(label), path("orthogroups_tradeoff.{pdf,png,svg}"), emit: tradeoff

    script:
    """
    # Evaluate BUSCO completeness and gene-occupancy threshold combinations.
    select_orthogroups.py \\
        --full-tables ${table.join(' ')} \\
        --output orthogroups \\
        --dataset-size "${dataset_size}"
    """

    stub:
    """
    command -v select_orthogroups.py >/dev/null

    touch orthogroups_busco-c_60_gene-occupancy_80.tsv
    touch orthogroups_summary.tsv

    touch orthogroups_heatmap.{pdf,png,svg}
    touch orthogroups_heatmap_zoom_50_100.{pdf,png,svg}
    touch orthogroups_tradeoff.{pdf,png,svg}
    """
}