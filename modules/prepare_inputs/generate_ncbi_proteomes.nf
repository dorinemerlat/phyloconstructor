process GENERATE_NCBI_PROTEOMES {
    tag "${specie}"
    // label 'agat'
    scratch false
    stageInMode 'copy'
    label 'retry_with_backoff'
    maxRetries 5
    cpus 4
    memory '32 GB'
    
    input:
    tuple val(taxid), val(specie), val(accession), path(genome), path(gff)

    output:
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}_canonical.gff"), emit: gff
    tuple val(taxid), val(specie), val(accession), path("${specie}_${accession}.fasta"), emit: fasta
    
    script:
    """
    # extract the longest isoform for each gene using AGAT
    module load agat
    agat_sp_keep_longest_isoform.pl --gff $gff -o ${specie}_${accession}_canonical.gff

    # extract proteins sequences
    agat_sp_extract_sequences.pl -g ${specie}_${accession}_canonical.gff -f $genome -o ${specie}_${accession}.fasta --protein
    """

    stub:
    """
    touch ${specie}_${accession}_canonical.gff ${specie}_${accession}.fasta
    """
}