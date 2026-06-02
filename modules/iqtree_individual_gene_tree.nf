process IQTREE_INDIVIDUAL_GENE_TREE {
    tag "${job_name}"
    cpus 20
    memory '20 GB'


    input:
    tuple val(job_name), val(outgroup), path(aln)

    output:
    tuple val(job_name), path("${outgroup}.treefile"), emit: treefile
    tuple val(job_name), path("${outgroup}.*") 

    script:
    """
    module load iqtree

    iqtree3 \\
        -s ${aln} \\
        -m MFP \\
        -B 1000 \\
        --alrt 1000 \\
        --bnni \\
        --prefix ${outgroup} \\
        -T ${task.cpus}
    """
}