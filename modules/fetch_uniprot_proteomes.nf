process FETCH_UNIPROT_PROTEOMES {
    tag ""

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.upids.tsv"), emit: ids 
    tuple val(taxid), path("${taxid}.upids.out"), emit: info 


    script:
    """
    curl -L -G "https://rest.uniprot.org/proteomes/stream" \
    --data-urlencode "compressed=false" \
    --data-urlencode "format=tsv" \
    --data-urlencode "fields=upid,organism,organism_id,protein_count,busco,cpd" \
    --data-urlencode "query=((taxonomy_id:${taxid}))" \
    -o ${taxid}.upids.out

    awk -F '\\t' '
        BEGIN { OFS="\\t" }
        NR == 1 {
            print "taxid","specie","sra"
            next
        }

        NR > 1 {
            specie_taxid = \$3
            specie_name = \$2
            upid = \$1

            sub(/\\(.*/, "", specie_name)
            gsub(/^ +| +\$/, "", specie_name)
            gsub(/ +/, "-", specie_name)
            specie_name = tolower(specie_name)

            print specie_taxid, specie_name, upid
        }
    ' ${taxid}.upids.out > ${taxid}.upids.tsv
    """
}