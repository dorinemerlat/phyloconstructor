process CLIPKIT {
    tag "${label}/${job_name}"
    scratch false

    input:
    tuple val(label), val(job_name), val(orthogroup), path(aln)

    output:
    tuple val(label), val(job_name), val(orthogroup), path("${orthogroup}_clean.aln")

    script:
    """
    module load clipkit
    clipkit ${aln} -m smart-gap -o ${orthogroup}_clean.aln
    """
}