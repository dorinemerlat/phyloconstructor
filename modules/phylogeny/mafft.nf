process IQTREE_SUPERMATRIX {
    tag "${label}/${job_name}"
    label 'iqtree'

    cpus 20
    memory { "${64 + (32 * (task.attempt - 1))} GB" }
    time '5d'

    input:
    tuple val(label), val(job_name), path(aln), path(partition)

    output:
    tuple val(label), val(job_name), path("${label}_${job_name}_supermatrix.treefile")
    tuple val(label), val(job_name), path("${label}_${job_name}_supermatrix.*")

    script:
    def prefix = "${label}_${job_name}_supermatrix"

    """
    # Infer a partitioned maximum-likelihood tree from the concatenated alignment.
    iqtree3 \\
        -s "${aln}" \\
        -p "${partition}" \\
        -m MFP+MERGE \\
        -B 1000 \\
        --alrt 1000 \\
        --bnni \\
        --prefix "${prefix}" \\
        -T ${task.cpus}
    """

    stub:
    def prefix = "${label}_${job_name}_supermatrix"

    """
    command -v iqtree3 >/dev/null

    touch "${prefix}.treefile"
    touch "${prefix}.log"
    """
}
