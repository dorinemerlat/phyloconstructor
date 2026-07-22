process SELECT_ORTHOGROUPS {
    tag "$label"
    // publishDir "${params.outdir}/orthogroups/$label",
    //     mode: 'copy',
    //     pattern: "{orthogroups_summary.tsv,orthogroups_heatmap.*,orthogroups_heatmap_zoom_50_100.*}"

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
    select_orthogroups.py \\
        --full-tables ${table.join(' ')} \\
        --output orthogroups \\
        --dataset-size ${dataset_size}
    """
}