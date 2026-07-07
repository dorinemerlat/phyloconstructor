process ASTRAL {
    tag "${label}/${job_name}"
    cpus 20
    memory { "${50 * task.attempt} GB" }
    scratch false

    input:
    tuple val(label), val(job_name), path(treefile, stageAs: "input/*")

    output:
    tuple val(label), val(job_name), path("${label}_${job_name}.astral.tree") 
    tuple val(label), val(job_name), path("${label}_${job_name}.astral.log") 

    script:
    """
    module load astral

    cat input/*.treefile > gene_trees.tree

    astral \\
        --input gene_trees.tree \\
        --output ${label}_${job_name}.astral.tree \\
        --branch-annotate 2 \
        > ${label}_${job_name}.astral.log 2>&1
    """
}