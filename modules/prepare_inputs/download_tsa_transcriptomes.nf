process DOWNLOAD_TSA_TRANSCRIPTOMES {
    tag "${specie}/${tsa}"
    label 'rna_tools'
    cache 'lenient'
    cpus 4
    memory '16 GB'
    time '12h'
    maxForks 2

    input:
    tuple val(taxid), val(specie), val(tsa)

    output:
    tuple val(taxid), val(specie), val(tsa), path("${specie}_${tsa}.fasta")

    script:
    """
    export TMPDIR="\$PWD/tmp_tsa_${tsa}"
    export TEMP="\$TMPDIR"
    export TMP="\$TMPDIR"

    mkdir -p "\$TMPDIR"

    # Download the TSA archive and convert it to FASTA.
    prefetch -f ALL "${tsa}" \\
        --output-directory .

    fasterq-dump \\
        --fasta \\
        --threads ${task.cpus} \\
        --mem 1000M \\
        --temp "\$TMPDIR" \\
        "${tsa}" \\
        -O . \\
        -o "${specie}_${tsa}.fasta"

    test -s "${specie}_${tsa}.fasta" || {
        echo "ERROR: TSA FASTA is empty for ${tsa}" >&2
        exit 1
    }

    rm -rf "${tsa}" "\$TMPDIR"
    """

    stub:
    """
    command -v prefetch >/dev/null
    command -v fasterq-dump >/dev/null

    touch "${specie}_${tsa}.fasta"
    """
}