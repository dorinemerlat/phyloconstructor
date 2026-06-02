process STRINGTIE {
    tag "${specie}/${data_id}"

    cpus 8
    memory { "${20 + (10 * (task.attempt - 1))} GB" }
    time '2d'
    maxRetries 5

    input:
    tuple val(taxid), val(specie), val(data_id), path(bam), path(genome)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_stringtie.gtf"), emit: gtf
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_stringtie_transcripts.fa"), emit: transcripts

    script:
    """
    module load stringtie/3.0.3 gffread/0.12.7

    stringtie \\
        ${bam} \\
        -p ${task.cpus} \\
        -o ${specie}_${data_id}_stringtie.gtf

    gffread \\
        ${specie}_${data_id}_stringtie.gtf \\
        -g ${genome} \\
        -w ${specie}_${data_id}_stringtie_transcripts.fa
    """
}