process ASTRAL {
    tag "${label}/${job_name}"
    label 'astral'

    cpus 1
    memory { "${50 * task.attempt} GB" }
    time '2d'

    input:
    tuple val(label), val(job_name), path(treefile, stageAs: "input/*")

    output:
    tuple val(label), val(job_name), path("${label}_${job_name}.astral.tree")
    tuple val(label), val(job_name), path("${label}_${job_name}.astral.log")

    script:
    """
    # Combine the individual gene trees into a single ASTRAL input file.
    cat input/*.treefile > gene_trees.tree

    # Infer the species tree and report local posterior probabilities.
    astral \\
        --input gene_trees.tree \\
        --output "${label}_${job_name}.astral.tree" \\
        --branch-annotate 2 \\
        > "${label}_${job_name}.astral.log" 2>&1
    """

    stub:
    """
    command -v astral >/dev/null

    touch "${label}_${job_name}.astral.tree"
    touch "${label}_${job_name}.astral.log"
    """
}