process CD_HIT_EST {
    tag "$id"
    cpus 5
    memory {"${800 + (200 * (task.attempt - 1))} MB"}
    time '1d'
    label 'phyloconstructor'
    maxRetries 5
    
    input:
    tuple val(id), val(input)

    output:
    tuple val(id), path("${id}_reformat.fasta")

    script:
    """
    memory=\$(echo ${task.memory.toMega()} | cut -d'.' -f1)
    cd-hit-est -i $input -o ${id}_reformat.fasta -c 0.95 -g 1 -T $task.cpus -M \${memory} 
    """

    stub:
    """
    touch ${id}_reformat.fasta
    """
}