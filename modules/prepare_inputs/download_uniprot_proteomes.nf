process DOWNLOAD_UNIPROT_PROTEOMES {
    tag "${specie}"

    input:
    tuple val(taxid), val(specie), val(upid)

    output:
    tuple val(taxid), val(specie), val(upid), path("${specie}_${upid}.faa")

    script:
    """
    curl --silent --fail --location \\
        --output ${specie}_${upid}.faa.gz \\
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3A${upid}%29"

    gunzip -c ${specie}_${upid}.faa.gz > ${specie}_${upid}.faa
    """
}