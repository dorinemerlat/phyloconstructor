process GROUPING_ORTHOLOGS_BY_SPECIES {
    tag "${taxid}"

    input:
    tuple val(taxid), path(busco_sequences), val(marker_number)

    output:
    tuple val(taxid), path("${taxid}/"), path("single_copy_orthologs.csv")

    script:
    """
    mkdir $taxid
    cp --update=none */*faa ${taxid}/. 2>/dev/null || true


    file_number=\$(ls ${taxid}/*.faa 2>/dev/null |wc -l || true)
    if [ "\$file_number" -gt 0 ]; then
        # rename all files finish by '.faa' to 'taxid.faa'
        for file in ${taxid}/*.faa; do
            mv "\$file" "\${file%.faa}_${taxid}.faa"
        done
    else 
        # if there is no faa file, create an empty file
        touch ${taxid}/0000_${taxid}.faa
    fi

    # calculate the percentage of single copy ortholog sequences
    echo "single_copy_number,single_copy_percentage,maker_number" > single_copy_orthologs.csv

    percentage=\$((\${file_number}*100/${marker_number}))
    echo "\${file_number},\${percentage},${marker_number}" >> single_copy_orthologs.csv
    """

    stub:
    """
    mkdir ${taxid}
    touch ${taxid}/0000.faa
    touch single_copy_orthologs.csv
    """
}