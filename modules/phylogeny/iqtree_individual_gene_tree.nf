process IQTREE_INDIVIDUAL_GENE_TREE {
    tag "${label}/${job_name}"
    cpus 10
    memory { "${20 * task.attempt} GB" }


    input:
    tuple val(label), val(job_name), val(orthogroup), path(aln)

    output:
    tuple val(label), val(job_name), path("${orthogroup}.treefile"), emit: treefile
    tuple val(label), val(job_name), path("${orthogroup}.*") 

    script:
    """
    module load iqtree

    iqtree3 \\
        -s ${aln} \\
        -m LG+G4 \\
        -B 1000 \\
        --bnni \\
        --prefix ${orthogroup} \\
        -T ${task.cpus}
    """
}