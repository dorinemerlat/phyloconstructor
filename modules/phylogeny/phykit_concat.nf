process PHYKIT_CONCAT {
    tag "${label}/${job_name}"
    scratch false
    stageInMode 'copy'

    memory { "${10 + (5 * (task.attempt - 1))} GB" }
    maxRetries 4
    errorStrategy { task.attempt <= 5 ? 'retry' : 'ignore' }

    input:
    tuple val(label), val(job_name), path(aln, stageAs: "input/*")

    output:
    tuple val(label), val(job_name), path("${job_name}.fa"), path("${job_name}.partition"), emit: main_output
    tuple val(label), val(job_name), path("${job_name}.occupancy"), emit: occupancy

    script:
    """
    module load phykit
    module load seqkit

    mkdir cleaned

    for f in input/*; do
        seqkit replace \\
            -p '^(.*?)_source=.*' \\
            -r '\$1' \\
            "\$f" > "cleaned/\$(basename "\$f")"
    done

    find cleaned -type f | sort > alignments.list

    phykit create_concat \\
        -a alignments.list \\
        -p ${job_name}

    seqkit seq -w 60 ${job_name}.fa > ${job_name}.fa.tmp
    mv ${job_name}.fa.tmp ${job_name}.fa

    test -s ${job_name}.fa || { echo "ERROR: ${job_name}.fa is empty" >&2; exit 1; }
    test -s ${job_name}.partition || { echo "ERROR: ${job_name}.partition is empty" >&2; exit 1; }
    test -s ${job_name}.occupancy || { echo "ERROR: ${job_name}.occupancy is empty" >&2; exit 1; }
    """
}
