process EXTRACT_BUSCO_FASTA {
    tag "${job_name}"
    scratch false

    input:
    tuple val(job_name), path(table) 
    //, path(busco_sequences, stageAs: "input/*")

    output:
    tuple val(job_name), path("orthogroups_${job_name}*.fasta")

    script:
    """
    extract_busco_fasta.py \\
        --table ${table} \\
        --output-prefix orthogroups_${job_name} \\
        --output-dir . \\
        --busco-sequences $projectDir/cache/busco/*/*/*/*fasta
    """
}