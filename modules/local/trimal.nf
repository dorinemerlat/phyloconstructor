process TRIMAL {
    tag "$id"
    cpus 20
    label 'trimal'
    
    input:
    tuple val(id), val(sco_name), path(aln) 

    output:
    tuple val(id), val(sco_name), path("${sco_name}_clean.aln")

    script:
    """
    trimal -in ${aln} -out ${sco_name}_clean.aln -automated1 -fasta -noallgaps
    """

    stub:
    """
    touch ${sco_name}.aln
    """
}