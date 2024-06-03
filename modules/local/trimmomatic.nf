process TRIMMOMATIC {
    tag "$id"
    cpus 40
    memory "20 GB"
    time '5d'
    label 'phyloconstructor'

    input:
    tuple val(id), path(fastq1), path(fastq2)

    output:
    tuple val(id), path("${id}_1_trim.fa"), path("${id}_2_trim.fa")

    script:
    """
    trimmomatic PE -threads 20 -basein ${fastq1} -baseout ${id}.trim \
        -baseout ${id}_paired.fq.gz ${id}_1_paired.fa ${id}_1_unpaired.fa ${id}_2_paired.fa ${id}_2_unpaired.fa \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    cat ${id}_1_paired.fa ${id}_1_unpaired.fa > ${id}_1_trim.fa
    cat ${id}_2_paired.fa ${id}_2_unpaired.fa > ${id}_2_trim.fa
    """

    stub:
    """
    touch ${id}_1_trim.fa ${id}_2_trim.fa
    """
}
