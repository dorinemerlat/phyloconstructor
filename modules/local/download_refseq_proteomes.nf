process DOWNLOAD_REFSEQ_PROTEOMES {
    tag "${taxid}"
    cache 'lenient'
    label 'phyloconstructor'

    input:
    val(taxid)

    output:
    tuple val(taxid), path{"${taxid}_refseq_prot.csv"}, path{"*.fna"}

    script:
    """
    datasets summary genome taxon $taxid > genome_summary.txt || true 

    echo "organism,refseq_prot_accession" > ${taxid}_refseq_prot.csv

    if [ -s genome_summary.txt ]; then
        datasets download genome taxon $taxid --include protein || true 
        
        unzip ncbi_dataset.zip

        # get all accession where proteins are available
        if [ -s ncbi_dataset/data/*/protein.faa ]; then
            for set in \$(ls ncbi_dataset/data/*/protein.faa); do
                name=\$(echo \$set |cut -d '/' -f 3)
                echo \$name >> ${taxid}_refseq_prot.accession
            done

            # get the taxid and accession of all genomes
            jq -r '[.organism.taxId, .accession] | @csv' ncbi_dataset/data/assembly_data_report.jsonl |sed 's/"//g' > ${taxid}_genomes.csv

            # find the organism name for each accession
            while read i; do
                organism=\$(grep \$i ${taxid}_genomes.csv| cut -d ',' -f 1)
                echo "\$organism,\${i}_prot" >> ${taxid}_refseq_prot.csv
            done < ${taxid}_refseq_prot.accession

            mv ncbi_dataset/data/\$name/protein.faa \${name}_prot.fna
        else 
            touch 0000_prot.fna
        fi

    else 
        touch 0000_prot.fna
    fi
    """

    stub:
    """
    touch ${taxid}_genomes.csv 0000.fna
    """
}