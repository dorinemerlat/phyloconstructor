process DOWNLOAD_NCBI_ASSEMBLIES {
    tag "${specie}/${accession}"
    label 'ncbi_datasets'
    cache 'lenient'
    memory '4 GB'
    time '4h'
    maxForks 20

    input:
    tuple val(taxid), val(specie), val(accession)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.fna"), path("${specie}_${accession}.gff")

    script:
    def download_dir = "${specie}_${accession}/ncbi_dataset/data/${accession}"

    """
    # Stagger requests to reduce simultaneous access to the NCBI API.
    sleep \$(shuf -i 1-20 -n 1)

    datasets download genome accession "${accession}" \\
        --include genome,gff3 \\
        --filename "${specie}_${accession}.zip"

    unzip -q "${specie}_${accession}.zip" \\
        -d "${specie}_${accession}"

    fna_file=\$(find "${download_dir}" -type f -name "*_genomic.fna" -print -quit)

    if [[ -z "\$fna_file" ]]
    then
        echo "ERROR: no genomic FASTA found for ${specie} ${accession}" >&2
        exit 1
    fi

    cp "\$fna_file" "${specie}_${accession}.fna"

    gff_file=\$(find "${download_dir}" -type f \\( -name "*.gff" -o -name "*.gff3" \\) -print -quit)

    if [[ -n "\$gff_file" ]]
    then
        cp "\$gff_file" "${specie}_${accession}.gff"
    else
        touch "${specie}_${accession}.gff"
    fi
    """

    stub:
    """
    command -v datasets >/dev/null
    command -v unzip >/dev/null

    touch "${specie}_${accession}.fna"
    touch "${specie}_${accession}.gff"
    """
}