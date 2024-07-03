process TRINITY {
    tag "$id"
    cpus 10
    memory {"${mem.toInteger() + 50 + (100 * (task.attempt -1))} GB"} // { "${memory.toInteger() + 20 + (10 * (task.attempt - 1))} GB" } // +10 GB for each retry
    time '5d'
    label 'trinity'
    maxRetries 5

    input:
    tuple val(id), path(fastq1), path(fastq2), val(mem)

    output:
    tuple val(id), path("trinity_${id}.Trinity.fasta")

    script:
    """
    replace_spaces_with_commas() {
        echo \$1 | sed 's/ /,/g'
    }

    fastq1_comma=\$(replace_spaces_with_commas $fastq1)
    fastq2_comma=\$(replace_spaces_with_commas $fastq2)

    # Add supplementary options to minimize RAM usage if the fastq files have multiple read sets
    supplementary_options=""
    if [[ \$fastq1_comma == *","* ]]; then
        supplementary_options="\${supplementary_options} --normalize_by_read_set"
    fi

    # if task.memory >= 50
    memory=\$(echo ${task.memory.toGiga()} | cut -d'.' -f1)
    if [[ \$memory -ge 50 ]]; then
        supplementary_options="\${supplementary_options} --no_parallel_norm_stats"
    fi

    Trinity --seqType fq --left \${fastq1_comma} --right \${fastq2_comma} --output trinity_${id} --CPU $task.cpus --max_memory \${memory}G --trimmomatic \$supplementary_options
    """

    stub:
    """
    touch trinity_${specie_name}.Trinity.fasta
    """
}
