process TRIMAL {
    tag "${job_name}"


    input:
    tuple val(job_name), val(orthogroup), path(aln) 

    output:
    tuple val(job_name), val(orthogroup), path("${orthogroup}_clean.aln") 

    script:
    """
    module load trimal
    trimal -in ${aln} -out ${orthogroup}_clean.aln -automated1 -fasta
    """
}