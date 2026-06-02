process DOWNLOAD_NCBI_ASSEMBLIES {
    tag "${specie}"
    cache 'lenient'
    maxForks 20

    input:
    tuple val(taxid), val(specie), val(accession)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.fna")

    script:
    """
    module load ncbi-datasets-cli

    sleep \$(shuf -i 1-20 -n 1)

    datasets download genome accession ${accession} \\
        --include genome \\
        --filename ${specie}_${accession}.genome.zip

    unzip -q ${specie}_${accession}.genome.zip -d ${specie}_${accession}.genome

    find ${specie}_${accession}.genome/ncbi_dataset/data/${accession} \\
        -name "*_genomic.fna" \\
        -exec cat {} \\; > ${specie}_${accession}.fna
    """
}