process DOWNLOAD_NCBI_PROTEOMES {
    tag "${specie}"
    cache 'lenient'
    maxForks 10

    input:
    tuple val(taxid), val(specie), val(accession)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.faa")

    script:
    """
    module load ncbi-datasets-cli

    sleep \$(shuf -i 1-20 -n 1)

    datasets download genome accession ${accession} \\
        --include protein \\
        --filename ${specie}_${accession}.protein.zip

    unzip -q ${specie}_${accession}.protein.zip -d ${specie}_${accession}.protein

    find ${specie}_${accession}.protein/ncbi_dataset/data/${accession} \\
        -name "protein.faa" \\
        -exec cat {} \\; > ${specie}_${accession}.faa
    """
}