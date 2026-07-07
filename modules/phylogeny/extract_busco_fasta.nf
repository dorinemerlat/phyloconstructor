process EXTRACT_BUSCO_FASTA {
    tag "${label}/${job_name}"
    scratch false

    input:
    tuple val(label), val(job_name), path(table), path(single_busco_sequences, stageAs: "single_sequences/*"), path(multi_busco_sequences, stageAs: "multi_sequences/*")

    output:
    tuple val(label), val(job_name), path("orthogroups_${job_name}*.fasta")

    script:
    """
    extract_busco_fasta.py \\
        --table ${table} \\
        --output-prefix orthogroups_${job_name} \\
        --output-dir . \\
        --busco-sequences single_sequences/* multi_sequences/*
    """
}