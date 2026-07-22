process FILTER_TRANSCRIPTS {
    tag "${specie}/${source}/${data_id}"
    label 'sequence_tools'

    cpus 20
    memory { "${100 * task.attempt} GB" }
    time '1d'
    maxRetries 5
    cache 'lenient'

    input:
    tuple val(taxid), val(specie), val(data_id), path(transcripts), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.min300.nr.fasta"), val(source)

    script:
    def prefix = "${specie}_${source}_${data_id}"

    """
    # Remove transcripts shorter than 300 nucleotides.
    seqkit seq \\
        -m 300 \\
        "${transcripts}" \\
        > "${prefix}.min300.fasta"

    # Cluster transcripts at 95% nucleotide identity.
    memory_mb=${task.memory.toMega() as long}

    cd-hit-est \\
        -i "${prefix}.min300.fasta" \\
        -o "${prefix}.min300.nr.fasta" \\
        -c 0.95 \\
        -g 1 \\
        -T ${task.cpus} \\
        -M "\$memory_mb"
    """

    stub:
    """
    command -v seqkit >/dev/null
    command -v cd-hit-est >/dev/null

    touch "${specie}_${source}_${data_id}.min300.nr.fasta"
    """
}