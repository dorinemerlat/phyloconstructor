process REFORMAT_FASTA {
    tag "${id}"

    input:
    tuple val(id), path(fasta) 

    output:
    tuple val(id), path("${id}_reformat.fasta")

    script:
    """
    sed -i "s|/|_|g" $fasta > ${id}_reformat.fasta
    """

    stub:
    """
    ls
    """
}