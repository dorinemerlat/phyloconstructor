process DOWNLOAD_TSA_TRANSCRIPTOMES {
    tag "${specie}"
    cache 'lenient'

    cpus 4
    memory '16 GB'
    maxForks 2

    input:
    tuple val(taxid), val(specie), val(tsa)

    output:
    tuple val(taxid), val(specie), val(tsa), path("${specie}_${tsa}.fasta")

    script:
    """
    module load sra-tools/3.4.1

    export TMPDIR=\$PWD/tmp_tsa_${tsa}
    export TEMP="\$TMPDIR"
    export TMP="\$TMPDIR"

    mkdir -p "\$TMPDIR"

    prefetch -f ALL ${tsa} --output-directory .

    fasterq-dump \\
        --fasta \\
        --threads ${task.cpus} \\
        --mem 1000M \\
        --temp "\$TMPDIR" \\
        ${tsa} \\
        -O . \\
        -o ${specie}_${tsa}.fasta

    rm -rf "\$TMPDIR"
    """

    stub:
    """
    touch ${specie}_${tsa}.fasta
    """
}