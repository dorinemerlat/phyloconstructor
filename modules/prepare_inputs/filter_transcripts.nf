process FILTER_TRANSCRIPTS {
    tag "${specie}/${source}"
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
    """
    module load seqkit cd-hit

    seqkit seq \
        -m 300 \
        "$transcripts" \
        > "${specie}_${source}_${data_id}.min300.fasta"
    
    memory=\$(echo ${task.memory.toMega()} | cut -d'.' -f1)
    cd-hit-est -i ${specie}_${source}_${data_id}.min300.fasta -o ${specie}_${source}_${data_id}.min300.nr.fasta -c 0.95 -g 1 -T ${task.cpus} -M \${memory} 
    """
}