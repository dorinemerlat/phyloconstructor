process GENERATE_NCBI_PROTEOMES {
    tag "${specie}/${accession}"
    label 'sequence_tools'
    cpus 4
    memory '32 GB'
    time '12h'
    maxRetries 5

    input:
    tuple val(taxid), val(specie), val(accession), path(genome), path(gff)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}_canonical.gff"), emit: gff
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.fasta"), emit: fasta

    script:
    """
    # Retain the longest transcript isoform for each annotated gene.
    agat_sp_keep_longest_isoform.pl \\
        --gff "${gff}" \\
        -o "${specie}_${accession}_canonical.gff"

    # Extract the corresponding protein sequences from the genome.
    agat_sp_extract_sequences.pl \\
        -g "${specie}_${accession}_canonical.gff" \\
        -f "${genome}" \\
        -o "${specie}_${accession}.fasta" \\
        --protein

    test -s "${specie}_${accession}.fasta" || {
        echo "ERROR: no protein sequence was generated for ${accession}" >&2
        exit 1
    }
    """

    stub:
    """
    command -v agat_sp_keep_longest_isoform.pl >/dev/null
    command -v agat_sp_extract_sequences.pl >/dev/null

    touch "${specie}_${accession}_canonical.gff"
    touch "${specie}_${accession}.fasta"
    """
}