process ASSEMBLE_SRA {
    tag "$taxid"
    cpus 40
    memory "20 GB"
    time '5d'
    label 'trinity'

    input:
    tuple val(taxid), path(fastq1), path(fastq2)

    output:
    tuple val(taxid), path("trinity_${taxid}.Trinity.fasta")

    script:
    """
    replace_spaces_with_commas() {
        echo \$1 | sed 's/ /,/g'
    }

    fastq1_comma=\$(replace_spaces_with_commas $fastq1)
    fastq2_comma=\$(replace_spaces_with_commas $fastq2)
    memory=\$(echo $task.memory | sed 's/ GB/G/')

    Trinity --seqType fq --left \${fastq1_comma} --right \${fastq2_comma} --output trinity_${taxid} --CPU $task.cpus --max_memory \$memory
    """

    stub:
    """
    touch trinity_${taxid}.Trinity.fasta
    """
}
