process TRIMAL {
    tag "${label}/${job_name}"
    scratch false
    
    input:
    tuple val(label), val(job_name), val(orthogroup), path(aln) 

    output:
    tuple val(label), val(job_name), val(orthogroup), path("${orthogroup}_clean.aln") 

    script:
    """
    module load trimal
    trimal -in ${aln} -out ${orthogroup}_clean.aln -gappyout -fasta
    """
}