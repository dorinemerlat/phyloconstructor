process DOWNLOAD_GENOMES {
    tag "${taxid}"
    cache 'lenient'
    label 'phyloconstructor'

    input:
    val(taxid)

    output:
    tuple val(taxid), path{"${taxid}_genomes.csv"}, path{"*.fna"}

    script:
    """
    datasets summary genome taxon $taxid > genome_summary.txt || true 

    if [ -s genome_summary.txt ]; then
        datasets download genome taxon $taxid --filename genome.zip
        unzip genome.zip

        echo "organism,genome_accession" > ${taxid}_genomes.csv
        jq -r '[.organism.taxId, .accession] | @csv' ncbi_dataset/data/assembly_data_report.jsonl |sed 's/"//g' >> ${taxid}_genomes.csv

        mv ncbi_dataset/data/*/*.fna .
    else 
        touch 0000.fna ${taxid}_genomes.csv
    fi
    """

    stub:
    """
    touch ${taxid}_genomes.csv 0000.fna
    """
}
