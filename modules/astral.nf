process ASTRAL {
    tag "${job_name}"
    cpus 20
    memory { "${50 * task.attempt} GB" }
    scratch false

    input:
    tuple val(job_name), path(treefile, stageAs: "input/*")

    output:
    tuple val(job_name), path("${job_name}.astral.tree") 
    tuple val(job_name), path("${job_name}.astral.log") 

    script:
    """
    module load astral

    cat input/*.treefile > gene_trees.tree

    astral \\
        --input gene_trees.tre \\
        --output ${job_name}.astral.tree \\
        --branch-annotate 2 \
        > ${job_name}.astral.log 2>&1
    """
}