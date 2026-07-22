process DOWNLOAD_NCBI_PROTEOMES {
    tag "${specie}/${accession}"
    label 'ncbi_datasets'
    cache 'lenient'
    memory '4 GB'
    time '4h'
    maxForks 10

    input:
    tuple val(taxid), val(specie), val(accession)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.faa")

    script:
    def download_dir = "${specie}_${accession}.protein/ncbi_dataset/data/${accession}"

    """
    # Stagger requests to reduce simultaneous access to the NCBI API.
    sleep \$(shuf -i 1-20 -n 1)

    datasets download genome accession "${accession}" \\
        --include protein \\
        --filename "${specie}_${accession}.protein.zip"

    unzip -q "${specie}_${accession}.protein.zip" \\
        -d "${specie}_${accession}.protein"

    find "${download_dir}" -type f -name "protein.faa" -exec cat {} + \\
        > "${specie}_${accession}.faa"

    if [[ ! -s "${specie}_${accession}.faa" ]]
    then
        echo "ERROR: no protein FASTA found for ${specie} ${accession}" >&2
        exit 1
    fi
    """

    stub:
    """
    command -v datasets >/dev/null
    command -v unzip >/dev/null

    touch "${specie}_${accession}.faa"
    """
}