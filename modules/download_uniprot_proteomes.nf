process DOWNLOAD_UNIPROT_PROTEOMES {
    tag "${specie}"

    input:
    tuple val(taxid), val(specie), val(upid)

    output:
    tuple val(taxid), val(specie), val(upid), path("${specie}_${upid}.faa")

    script:
    """
    curl --silent --output ${specie}_${upid}.faa "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3A${upid}%29%29"
    """
}