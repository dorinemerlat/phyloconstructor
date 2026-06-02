process CD_HIT_EST {
    tag "${specie}/${source}/${data_id}"
    cpus 15
    memory { "${150 + (20 * (task.attempt - 1))} GB" }
    time '1d'
    maxRetries 5
    
    input:
   tuple val(taxid), val(specie), val(data_id), path(input), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.nr.fa"), val(source)

    script:
    """
    module load cd-hit
    memory=\$(echo ${task.memory.toMega()} | cut -d'.' -f1)
    cd-hit-est -i $input -o ${specie}_${source}_${data_id}.nr.fa -c 0.95 -g 1 -T $task.cpus -M \${memory} 
    """
}