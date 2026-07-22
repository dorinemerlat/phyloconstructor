process TRANSDECODER {
    tag "${specie}/${source}/${data_id}"
    label 'transdecoder'
    cpus 2
    memory { "${20 + (4 * (task.attempt - 1))} GB" }
    time '20h'
    maxRetries 5

    input:
    tuple val(taxid), val(specie), val(data_id), path(input), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.transdecoder.pep"), val(source)

    script:
    """
    # Identify candidate long open reading frames.
    TransDecoder.LongOrfs \\
        -t "${input}"

    # Retain the best predicted coding sequence per transcript.
    TransDecoder.Predict \\
        -t "${input}" \\
        --single_best_only \\
        --no_refine_starts

    mv "${input}.transdecoder.pep" \\
        "${specie}_${source}_${data_id}.transdecoder.pep"
    """

    stub:
    """
    command -v TransDecoder.LongOrfs >/dev/null
    command -v TransDecoder.Predict >/dev/null

    touch "${specie}_${source}_${data_id}.transdecoder.pep"
    """
}