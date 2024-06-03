process DOWNLOAD_PROTEOMES {
    tag "${taxid}"
    cache 'lenient'
    label 'jq'
    
    input:
    val(taxid)

    output:
    tuple val(taxid), path("${taxid}_proteomes.csv"), path("*.faa")

    script:
    """
    # get all proteome IDs in the given clade
    curl --silent --output ${taxid}_UPIDs.list "https://rest.uniprot.org/proteomes/stream?format=list&query=%28%28taxonomy_id%3A${taxid}%29%29"

    echo "organism,genome_accession" > ${taxid}_proteomes.csv

    if [ -s ${taxid}_UPIDs.list ]; then
        cat ${taxid}_UPIDs.list | while read upid; do
            # get all proteins in a proteome
            curl --silent --output \${upid}.faa "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A\${upid}%29%29"

            # get information about the proteome
            curl --silent --output \${upid}_info.json "https://www.ebi.ac.uk/proteins/api/proteomes/\${upid}"

            jq -r '[.taxonomy, .upid] | @csv' \${upid}_info.json | sed 's/"//g' >> ${taxid}_proteomes.csv
        done
    else 
        touch ${taxid}_proteomes.csv 0000.faa
    fi
    """

    stub:
    """
    touch ${taxid}_proteomes.csv 0000.faa
    """
}
