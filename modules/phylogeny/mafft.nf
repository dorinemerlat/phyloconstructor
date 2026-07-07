process MAFFT {
    tag "${label}/${job_name}"
    cpus 4
    memory { "${15 + (20 * (task.attempt - 1))} GB" }

    input:
    tuple val(label), val(job_name), val(orthogroup), path(fasta) 
    //, path(busco_sequences, stageAs: "input/*")

    output:
    tuple val(label), val(job_name), val(orthogroup), path("${orthogroup}.aln") 

    script:
    """
    module load mafft
    mafft --localpair --maxiterate 1000 --reorder --thread ${task.cpus} ${fasta} > ${orthogroup}.aln
    """
}