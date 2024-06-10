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
    TrimmomaticPE -threads ${task.cpus} -basein ${fastq1} ${fastq2} \
        -baseout ${id}_1P.fq.gz ${id}_1P.fa.gz ${id}_1U.fa.gz ${id}_2P.fq.gz ${id}_2U.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    zcat ${id}_1P.fa.gz ${id}_1U.fa.gz > ${id}_1_trim.fa
    zcat ${id}_2P.fq.gz ${id}_2U.fq.gz > ${id}_2_trim.fa
    """

    stub:
    """
    touch ${id}_1_trim.fa ${id}_2_trim.fa
    """
}
