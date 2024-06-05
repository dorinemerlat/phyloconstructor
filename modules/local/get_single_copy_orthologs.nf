process GET_SINGLE_COPY_ORTHOLOGS {
    tag "${taxid}"
    scratch false
    
    input:
    tuple val(taxid), path(busco_sequences), val(marker_number)

    output:
    tuple val(taxid), path("single_copy_orthologs_${taxid}/*.faa"), path("single_copy_orthologs.csv")

    script:
    """
    mkdir -p single_copy_orthologs_${taxid}

    for busco_res in ${busco_sequences}; do
        for sequence in \$(ls \${busco_res}); do

            if [ ! -f single_copy_orthologs/\${sequence} ]; then
                name=\$(basename -s .tar.gz "\$sequence")
                cp \${busco_res}/\${sequence} single_copy_orthologs_${taxid}
            fi
        done
    done

    file_number=\$(ls single_copy_orthologs_${taxid}/*faa |wc -l)
    echo "single_copy_number,single_copy_percentage,maker_number" > single_copy_orthologs.csv
    percentage=\$((\${file_number}*100/${marker_number}))
    echo "\${file_number},\${percentage},${marker_number}" >> single_copy_orthologs.csv
    """

    stub:
    """
    mkdir -p ${taxid}
    touch ${taxid}/0000.faa
    """
}