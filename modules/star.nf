process STAR {
    tag "${specie}/${data_id}"

    cpus 20
    time '5d'
    label 'varus'
    scratch false
    memory { (100 * task.attempt).GB }
    stageInMode 'copy'
    maxRetries 10

    input:
    tuple val(taxid), val(specie), val(data_id), path(reads1), path(reads2), path(genome)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_star.bam"), emit: bam
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_star_Log.out"), emit: star_log
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_star_Log.progress.out"), emit: star_progress
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_star_Log.final.out"), emit: star_log_final
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${data_id}_star_mapped_read_count.txt"), emit: mapped_count

    script:
    def r1_csv = reads1 instanceof List ? reads1.sort { it.name }.join(',') : reads1.toString()
    def r2_csv = reads2 instanceof List ? reads2.sort { it.name }.join(',') : reads2.toString()

    """
    mkdir -p genome

    STAR \\
        --runThreadN 4 \\
        --runMode genomeGenerate \\
        --outTmpDir STARtmp_index \\
        --genomeDir genome \\
        --genomeFastaFiles ${genome} \\
        --limitGenomeGenerateRAM ${task.memory.toBytes()}

    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir genome \\
        --readFilesIn ${r1_csv} ${r2_csv} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${specie}_${data_id}_star_ \\
        --outSAMstrandField intronMotif \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM ${task.memory.toBytes()}

    mv ${specie}_${data_id}_star_Aligned.sortedByCoord.out.bam ${specie}_${data_id}_star.bam

    samtools view -c -F 4 ${specie}_${data_id}_star.bam > ${specie}_${data_id}_star_mapped_read_count.txt
    """
}