process DOWNLOAD_SRA {
    tag "$taxid"
    cache 'lenient'
    label 'phyloconstructor'

    input:
    tuple val(taxid), val(ncbi_api_key) 

    output:
    tuple val(taxid), path("${taxid}_sra.csv"), path("*_1.fastq"), path("*_2.fastq")

    script:
    """
    export NCBI_API_KEY=$ncbi_api_key

    { # try a first time
         esearch -db sra -query '((((txid${taxid}[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \
            > esearch.out

    } || { # try a second time
        sleep \$(shuf -i 5-30 -n 1)
        esearch -db sra -query '((((txid${taxid}[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \
            > esearch.out
    }

    count=\$(grep "<Count>" esearch.out |cut -d '>' -f 2 |cut -d '<' -f1)
    if [[ \$count != 0 ]] ; then
        efetch -format runinfo < esearch.out | awk -F','  'BEGIN {OFS=","} {print \$28,\$1}' > ${taxid}_sra.csv

        tail -n +2 ${taxid}_sra.csv | cut -d ',' -f 2 > sra.accession
        prefetch -f ALL --option-file sra.accession

        srr_directories=\$(find . -type d -name 'SRR*')

        for srr in \${srr_directories}; do
            fasterq-dump --split-files \${srr}
        done
    else
        touch 0000_1.fastq 0000_2.fastq
    fi
    """

    stub:
    """
    touch ${taxid}_sra.csv 0000_1.fastq 0000_2.fastq
    """
}
