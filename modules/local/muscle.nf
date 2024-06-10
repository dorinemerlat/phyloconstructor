process MUSCLE {
    tag "$id"
    cpus 20
    label 'phyloconstructor'

    input:
    tuple val(id), val(sco_name), path(input) 

    output:
    tuple val(id), val(sco_name), path("${sco_name}.aln")

    script:
    """
    muscle -super5 $input -output ${sco_name}.aln  -threads ${task.cpus}
    """

    stub:
    """
    touch ${sco_name}.aln
    """
}