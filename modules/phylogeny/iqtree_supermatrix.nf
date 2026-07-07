process IQTREE_SUPERMATRIX {
    tag "${label}/${job_name}"
    cpus 20
    memory { "${50 * task.attempt} GB" }
    scratch false
    time '2d'

    input:
    tuple val(label), val(job_name), path(aln), path(partition)

    output:
    tuple val(label), val(job_name), path("${label}_${job_name}_supermatrix.treefile") 
    tuple val(label), val(job_name), path("${label}_${job_name}_supermatrix.*") 

    script:
    """
    module load iqtree

    iqtree3 \\
        -s ${aln} \\
        -p ${partition} \\
        -m MFP+MERGE \\
        -B 1000 \\
        --alrt 1000 \\
        --bnni \\
        --prefix ${label}_${job_name}_supermatrix \\
        -T ${task.cpus} 
    """
}