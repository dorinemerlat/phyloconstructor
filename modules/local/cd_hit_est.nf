process CD_HIT_EST {
    tag "$id"
    cpus 20
    time '1d'
    label 'phyloconstructor'

    input:
    tuple val(id), val(input)

    output:
    tuple val(id), path("${id}_reformat.fasta")

    script:
    """
    cd-hit-est -i $input -o ${id}_reformat.fasta -c 0.95 -g 1 -T $task.cpus 
    """

    stub:
    """
    touch ${id}_reformat.fasta
    """
}