process TRINITY {
    tag "$id"
    cpus 40
    memory "50 GB"
    time '5d'
    label 'trinity'
    maxRetries 0

    input:
    tuple val(id), path(fastq1), path(fastq2)

    output:
    tuple val(id), path("trinity_${id}.Trinity.fasta")

    script:
    """
    replace_spaces_with_commas() {
        echo \$1 | sed 's/ /,/g'
    }

    fastq1_comma=\$(replace_spaces_with_commas $fastq1)
    fastq2_comma=\$(replace_spaces_with_commas $fastq2)
    memory='${task.memory}'
    memory=\${memory%B}
    memory=\${memory// /}

    Trinity --seqType fq --left \${fastq1_comma} --right \${fastq2_comma} --output trinity_${id} --CPU $task.cpus --max_memory \$memory --trimmomatic
    """

    stub:
    """
    touch trinity_${specie_name}.Trinity.fasta
    """
}
