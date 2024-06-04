process GROUP_BY_ORTHOLOGS {
    tag "${id}"

    input:
    tuple val(id), path(sequences) 

    output:
    tuple val(id), path('single_copy_orthologs/')

    script:
    """
    mkdir single_copy_orthologs

    for directory in $(ls sequences); do
        for sequence in \$(ls \${directory}); do
            name=\$(basename -s .tar.gz "$sequence")

                mkdir -p single_copy_orthologs/$name
                cp sequences/$directory/$sequence single_copy_orthologs/${name}_${sequence}.faa
        
        done
    done
    """

    stub:
    """
    mdkir single_copy_orthologs
    touch single_copy_orthologs/0000.faa
    """
}