process FETCH_TSA_TRANSCRIPTOMES {
    tag "${taxid}"
    label 'entrez_direct'
    memory '2 GB'
    time '2h'

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.tsa_ids.tsv"), emit: ids
    tuple val(taxid), path("${taxid}.tsa_ids.out"), emit: info

    script:
    """
    esearch -db nuccore \\
        -query "(txid${taxid}[Organism:exp]) AND \\"tsa master\\"[Properties]" \\
        > esearch.out

    count=\$(grep -o '<Count>[0-9]*</Count>' esearch.out | head -n 1 | sed 's/<[^>]*>//g')
    count=\${count:-0}

    printf "taxid\\tspecie\\ttsa\\n" \\
        > "${taxid}.tsa_ids.tsv"

    if [[ "\$count" != "0" ]]
    then
        efetch -format gb < esearch.out \\
            > "${taxid}.tsa_ids.out"

        grep "^ACCESSION" "${taxid}.tsa_ids.out" \\
            | awk '{ print \$2 }' \\
            > "${taxid}.tsa.list"

        grep '/db_xref="taxon:' "${taxid}.tsa_ids.out" \\
            | sed 's/.*taxon://' \\
            | cut -d '"' -f 1 \\
            > "${taxid}.organism_taxid.list"

        grep "^  ORGANISM" "${taxid}.tsa_ids.out" \\
            | sed 's/^  ORGANISM  //' \\
            | awk '
                {
                    specie = \$0
                    sub(/\\(.*/, "", specie)
                    gsub(/^ +| +\$/, "", specie)
                    gsub(/ +/, "-", specie)
                    print tolower(specie)
                }
            ' \\
            > "${taxid}.specie.list"

        paste \\
            "${taxid}.organism_taxid.list" \\
            "${taxid}.specie.list" \\
            "${taxid}.tsa.list" \\
            >> "${taxid}.tsa_ids.tsv"
    else
        touch "${taxid}.tsa_ids.out"
    fi
    """

    stub:
    """
    command -v esearch >/dev/null
    command -v efetch >/dev/null

    touch "${taxid}.tsa_ids.tsv"
    touch "${taxid}.tsa_ids.out"
    """
}