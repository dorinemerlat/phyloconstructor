process AMAS {
    tag "${job_name}"
    scratch false

    input:
    tuple val(job_name), path(aln, stageAs: "input/*") 

    output:
    tuple val(job_name), path("${job_name}_supermatrix.aln"), path("${job_name}_partition.txt")

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