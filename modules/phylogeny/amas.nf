process AMAS {
    tag "${label}/${job_name}"
    scratch false
    memory { "${30 + (10 * (task.attempt - 1))} GB" }
    maxRetries = 4
    errorStrategy = { task.attempt <= 5 ? 'retry' : 'ignore' }

    input:
    tuple val(label), val(job_name), path(aln, stageAs: "input/*") 

    output:
    tuple val(label), val(job_name), path("${job_name}_supermatrix.aln"), path("${job_name}_partition.txt")

    script:
    """
    python /shared/projects/metainvert/phyloconstructor2/AMAS/amas/AMAS.py concat \\
        -i input/* \\
        -f fasta \\
        -d aa \\
        -t ${job_name}_supermatrix.aln \
        -p ${job_name}_partition.txt
    """
}