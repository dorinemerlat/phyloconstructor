process FETCH_TSA_TRANSCRIPTOMES {
    tag ""

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.tsa_ids.tsv"), emit: ids
    tuple val(taxid), path("${taxid}.tsa_ids.out"), emit: info

    script:
    """
    module load entrez-direct/22.4

    esearch -db nuccore \\
        -query "(txid${taxid}[Organism:exp]) AND \\"tsa master\\"[Properties]" \\
        > esearch.out

    count=\$(grep "<Count>" esearch.out | cut -d '>' -f 2 | cut -d '<' -f 1)

    echo -e "taxid\\tspecie\\ttsa" > ${taxid}.tsa_ids.tsv

    if [[ "\$count" != "0" ]]; then
        efetch -format gb < esearch.out > ${taxid}.tsa_ids.out

        grep "^ACCESSION" ${taxid}.tsa_ids.out \\
            | awk '{ print \$2 }' \\
            > ${taxid}.tsa.list

        grep "/db_xref=\\"taxon:" ${taxid}.tsa_ids.out \\
            | sed 's/.*taxon://' \\
            | cut -d '"' -f 1 \\
            > ${taxid}.organism_taxid.list

        grep "^  ORGANISM" ${taxid}.tsa_ids.out \\
            | sed 's/^  ORGANISM  //' \\
            | awk '
                BEGIN { OFS="\\t" }
                {
                    specie = \$0
                    sub(/\\(.*/, "", specie)
                    gsub(/^ +| +\$/, "", specie)
                    gsub(/ +/, "-", specie)
                    specie = tolower(specie)
                    print specie
                }
            ' > ${taxid}.specie.list

        paste \\
            ${taxid}.organism_taxid.list \\
            ${taxid}.specie.list \\
            ${taxid}.tsa.list \\
            >> ${taxid}.tsa_ids.tsv
    else
        echo "taxid,specie,tsa" >  ${taxid}.tsa_ids.out
    fi
    """
}