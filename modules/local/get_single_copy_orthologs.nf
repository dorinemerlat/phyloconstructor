process GET_SINGLE_COPY_ORTHOLOGS {
    tag "${taxid}"

    input:
    tuple val(taxid), va(id), path(busco_sequences), val(marker_number)

    output:
    tuple val(taxid), path("taxid/*.faa"), path("single_copy_orthologs.csv")

    script:
    """
    mkdir -p ${taxid}

    for busco_res in ${busco_sequences}); do
        for sequence in \$(ls \${busco_res}/single_copy_busco_sequences); do

            if [ ! -f single_copy_orthologs/\${sequence} ]; then
                name=\$(basename -s .tar.gz "$sequence")
                cp \${busco_res}/single_copy_busco_sequences/\${sequence} single_copy_orthologs
            if

        done
    done

    file_number=$(ls $taxid/*faa)
    echo "single_copy_number,single_copy_percentage,maker_number" > single_copy_orthologs.csv
    percentage=\$((\${file_number}*100/${marker_number}))
    echo "${file_number},${percentage},${marker_number}" >> single_copy_orthologs.csv
    """

    stub:
    """
    mkdir -p ${taxid}
    touch ${taxid}/0000.faa
    """
}