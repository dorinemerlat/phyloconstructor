process TRIMAL {
    tag "${label}/${job_name}/${orthogroup}"
    label 'trimal'
    memory '2 GB'
    time '2h'

    input:
    tuple val(label), val(job_name), val(orthogroup), path(aln)

    output:
    tuple val(label), val(job_name), val(orthogroup), path("${orthogroup}_clean.aln")

    script:
    """
    # Remove poorly aligned positions using trimAl's gappyout heuristic.
    trimal \\
        -in "${aln}" \\
        -out "${orthogroup}_clean.aln" \\
        -gappyout \\
        -fasta
    """

    stub:
    """
    command -v trimal >/dev/null

    touch "${orthogroup}_clean.aln"
    """
}