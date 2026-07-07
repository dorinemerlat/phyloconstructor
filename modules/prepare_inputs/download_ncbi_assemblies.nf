process DOWNLOAD_NCBI_ASSEMBLIES {
    tag "${specie}"
    cache 'lenient'
    maxForks 20
    scratch false

    input:
    tuple val(taxid), val(specie), val(accession)

    output:
    tuple val(taxid), val(specie), val(accession),
          path("${specie}_${accession}.fna"),
          path("${specie}_${accession}.gff")

    script:
    download_dir = "${specie}_${accession}/ncbi_dataset/data/${accession}"

    """
    module load ncbi-datasets-cli

    sleep \$(shuf -i 1-20 -n 1)

    datasets download genome accession ${accession} \\
        --include genome,gff3 \\
        --filename ${specie}_${accession}.zip

    unzip -q ${specie}_${accession}.zip -d ${specie}_${accession}

    fna_file=\$(find ${download_dir} -type f -name "*_genomic.fna" | head -n 1)

    if [ -z "\${fna_file}" ]; then
        echo "ERROR: No genomic FASTA found for ${specie} ${accession}" >&2
        exit 1
    fi

    cp "\${fna_file}" ${specie}_${accession}.fna

    gff_file=\$(find ${download_dir} -type f \\( -name "*.gff" -o -name "*.gff3" \\) | head -n 1)

    if [ -n "\${gff_file}" ]; then
        cp "\${gff_file}" ${specie}_${accession}.gff
    else
        touch ${specie}_${accession}.gff
    fi
    """
}