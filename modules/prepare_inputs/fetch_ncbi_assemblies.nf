process FETCH_NCBI_ASSEMBLIES {
    tag "${taxid}"
    label 'ncbi_datasets'
    memory '4 GB'
    time '4h'

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.ncbi_accessions.tsv"), emit: assembly_ids
    tuple val(taxid), path("${taxid}.ncbi_accessions_for_proteins.tsv"), emit: proteomes_ids
    tuple val(taxid), path("${taxid}.ncbi_accessions.out"), emit: info

    script:
    """
    # Retrieve the complete NCBI assembly summary.
    datasets summary genome taxon "${taxid}" --as-json-lines \\
        | dataformat tsv genome \\
        > "${taxid}.ncbi_accessions.out"

    # Retrieve the fields required to construct the pipeline input tables.
    datasets summary genome taxon "${taxid}" --as-json-lines \\
        | dataformat tsv genome \\
            --fields organism-tax-id,organism-name,accession,annotinfo-featcount-gene-protein-coding \\
        > "${taxid}.ncbi_accessions.tsv.tmp"

    awk -F '\\t' '
        BEGIN { OFS="\\t" }

        NR == 1 {
            print "taxid","specie","genome"
            next
        }

        {
            specie_taxid = \$1
            specie_name  = \$2
            accession    = \$3

            sub(/\\(.*/, "", specie_name)
            gsub(/^ +| +\$/, "", specie_name)
            gsub(/ +/, "-", specie_name)
            specie_name = tolower(specie_name)

            print specie_taxid, specie_name, accession
        }
    ' "${taxid}.ncbi_accessions.tsv.tmp" \\
        > "${taxid}.ncbi_accessions.tsv"

    awk -F '\\t' '
        BEGIN { OFS="\\t" }

        NR == 1 {
            print "taxid","specie","genome"
            next
        }

        {
            protein_coding_genes = \$4

            if (protein_coding_genes == "" || protein_coding_genes == "0")
                next

            specie_taxid = \$1
            specie_name  = \$2
            accession    = \$3

            sub(/\\(.*/, "", specie_name)
            gsub(/^ +| +\$/, "", specie_name)
            gsub(/ +/, "-", specie_name)
            specie_name = tolower(specie_name)

            print specie_taxid, specie_name, accession
        }
    ' "${taxid}.ncbi_accessions.tsv.tmp" \\
        > "${taxid}.ncbi_accessions_for_proteins.tsv"
    """

    stub:
    """
    command -v datasets >/dev/null
    command -v dataformat >/dev/null

    touch "${taxid}.ncbi_accessions.tsv"
    touch "${taxid}.ncbi_accessions_for_proteins.tsv"
    touch "${taxid}.ncbi_accessions.out"
    """
}