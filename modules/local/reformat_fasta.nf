process REFORMAT_FASTA {
    tag "${id}"
    label 'bioawk'

    input:
    tuple val(id), path(fasta) 

    output:
    tuple val(id), path("${id}_reformat.fasta")

    script:
    """
        bioawk -c fastx '\$name !~ /isoform/ {print ">"\$name; print \$seq}' $fasta | fold -w 60 \
            | sed "s|/|_|g"  \
            > ${id}_reformat.fasta
    """

    stub:
    """
    touch ${id}_reformat.fasta
    """
}