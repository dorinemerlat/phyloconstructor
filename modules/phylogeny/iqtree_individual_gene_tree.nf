process IQTREE_INDIVIDUAL_GENE_TREE {
    tag "${label}/${job_name}/${orthogroup}"
    label 'iqtree'
    cpus 10
    memory { "${20 * task.attempt} GB" }
    time '1d'

    input:
    tuple val(label), val(job_name), val(orthogroup), path(aln)

    output:
    tuple val(label), val(job_name), path("${orthogroup}.treefile"), emit: treefile
    tuple val(label), val(job_name), path("${orthogroup}.*")

    script:
    """
    # Infer a maximum-likelihood tree and estimate ultrafast bootstrap support.
    iqtree3 \\
        -s "${aln}" \\
        -m MFP+MERGE \\
        -B 1000 \\
        --bnni \\
        --prefix "${orthogroup}" \\
        -T ${task.cpus}
    """

    stub:
    """
    command -v iqtree3 >/dev/null

    touch "${orthogroup}.treefile"
    touch "${orthogroup}.log"
    """
}