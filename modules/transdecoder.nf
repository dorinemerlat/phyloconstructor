process TRANSDECODER {
    tag "${specie}"
    cpus 2
    memory { "${4 + (4 * (task.attempt - 1))} GB" }
    time '20h'
    maxRetries 5

    input:
    tuple val(taxid), val(specie), val(data_id), path(input), val(source)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.transdecoder.pep"), val(source)

    script:
    """
    module load transdecoder/5.7.0

    TransDecoder.LongOrfs -t ${input}
    TransDecoder.Predict -t ${input} --single_best_only --no_refine_starts

    mv ${input}.transdecoder.pep ${specie}_${source}_${data_id}.transdecoder.pep
    """
}