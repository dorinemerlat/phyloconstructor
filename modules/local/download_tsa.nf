process DOWNLOAD_TSA {
    tag "${taxid}"
    cache 'lenient'
    label 'ncbi'

    input:
    tuple val(taxid), val(ncbi_api_key) 

    output:
    tuple val(taxid), path("${taxid}_tsa.csv"), path("*.fasta")

    script:
    """
    export NCBI_API_KEY=$ncbi_api_key

    { # try a first time
        esearch -db nuccore -query '(txid${taxid}[Organism:exp]) AND "tsa master"[Properties]' > esearch.out

    } || { # try a second time
        sleep \$(shuf -i 5-30 -n 1)
        esearch -db nuccore -query '(txid${taxid}[Organism:exp]) AND "tsa master"[Properties]' > esearch.out
    }

    echo "organism,tsa" > ${taxid}_tsa.csv
    count=\$(grep "<Count>" esearch.out |cut -d '>' -f 2 |cut -d '<' -f1)
    if [[ \$count != 0 ]] ; then
        efetch -format gb < esearch.out > efetch.out
        grep "^ACCESSION" efetch.out |awk '{ print \$2}' > ${taxid}_tsa.list
        grep -o "/db_xref=" efetch.out | cut -d':' -f2 | cut -d'\"' -f1 > ${taxid}_organism.list

        paste -d ',' ${taxid}_organism.list ${taxid}_tsa.list >> ${taxid}_tsa.csv

        tail -n +2 ${taxid}_tsa.csv | cut -d ',' -f 2 > tsa.accession
        prefetch --option-file tsa.accession

        srr_directories=\$(find . -type f -regex "\\./[A-Z0-9]*")
        for srr in \${srr_directories} ; do
            fasterq-dump --fasta \${srr}
        done
    
    else
        touch 0000.fasta
    fi
    """

    stub:
    """
    touch 0000.fasta
    """
}