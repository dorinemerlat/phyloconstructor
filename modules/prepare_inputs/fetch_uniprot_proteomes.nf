process FETCH_UNIPROT_PROTEOMES {
    tag "${taxid}"
    memory '2 GB'
    time '2h'

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.upids.tsv"), emit: ids
    tuple val(taxid), path("${taxid}.upids.out"), emit: info

    script:
    """
    # Retrieve the reference proteomes associated with the requested taxon.
    curl \\
        --silent \\
        --show-error \\
        --fail \\
        --location \\
        --get \\
        --retry 5 \\
        --retry-all-errors \\
        "https://rest.uniprot.org/proteomes/stream" \\
        --data-urlencode "compressed=false" \\
        --data-urlencode "format=tsv" \\
        --data-urlencode "fields=upid,organism,organism_id,protein_count,busco,cpd" \\
        --data-urlencode "query=((taxonomy_id:${taxid}) AND (proteome_type:REFERENCE))" \\
        --output "${taxid}.upids.out"

    if [[ ! -s "${taxid}.upids.out" ]]
    then
        echo "ERROR: UniProt returned an empty file for taxid ${taxid}" >&2
        exit 1
    fi

    awk -F '\\t' '
        BEGIN { OFS="\\t" }

        NR == 1 {
            print "taxid","specie","upid"
            next
        }

        {
            specie_taxid = \$3
            specie_name  = \$2
            upid         = \$1

            sub(/\\(.*/, "", specie_name)
            gsub(/^ +| +\$/, "", specie_name)
            gsub(/ +/, "-", specie_name)
            specie_name = tolower(specie_name)

            print specie_taxid, specie_name, upid
        }
    ' "${taxid}.upids.out" \\
        > "${taxid}.upids.tsv"
    """

    stub:
    """
    command -v curl >/dev/null

    touch "${taxid}.upids.tsv"
    touch "${taxid}.upids.out"
    """
}