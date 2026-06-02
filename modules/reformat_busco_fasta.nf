process REFORMAT_BUSCO_FASTA {
    tag "${busco_threshold}/${gene_occupancy_threshold}"
    scratch false

    input:
    tuple val(busco_threshold), val(gene_occupancy_threshold), path(table), path(fasta_dir, stageAs: "input/*")

    output:
    tuple val(busco_threshold), val(gene_occupancy_threshold), path("orthogroups_busco-c_${busco_threshold}_gene-occupancy_${gene_occupancy_threshold}*.fasta")

    script:
    """
    extract_busco_fasta.py \\
        --table ${table} \\
        --output-prefix orthogroups_busco-c_${busco_threshold}_gene-occupancy_${gene_occupancy_threshold} \\
        --busco-sequences input/*
    """
}