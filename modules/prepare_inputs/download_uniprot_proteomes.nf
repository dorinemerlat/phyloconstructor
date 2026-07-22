process DOWNLOAD_UNIPROT_PROTEOMES {
    tag "${specie}/${upid}"
    cache 'lenient'
    memory '2 GB'
    time '4h'
    maxForks 10

    input:
    tuple val(taxid), val(specie), val(upid)

    output:
    tuple val(taxid), val(specie), val(upid), path("${specie}_${upid}.faa")

    script:
    """
    # Download the compressed UniProt proteome.
    curl \\
        --silent \\
        --show-error \\
        --fail \\
        --location \\
        --retry 5 \\
        --retry-all-errors \\
        --output "${specie}_${upid}.faa.gz" \\
        "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3A${upid}%29"

    gzip -cd "${specie}_${upid}.faa.gz" \\
        > "${specie}_${upid}.faa"

    test -s "${specie}_${upid}.faa" || {
        echo "ERROR: UniProt returned an empty proteome for ${upid}" >&2
        exit 1
    }
    """

    stub:
    """
    command -v curl >/dev/null
    command -v gzip >/dev/null

    touch "${specie}_${upid}.faa"
    """
}