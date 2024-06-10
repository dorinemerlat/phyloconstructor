process GROUPING_ORTHOLOGS_BY_SCO_GROUPS {
    tag "$specie_threshold"

    input:
    tuple val(specie_threshold), path(sequences_dir) 

    output:
    tuple val(specie_threshold), path("sco_*.faa")

    script:
    """
    for sequence in */*faa; do
        name=\$(basename -s ".faa" \$sequence)
        sco=\$(echo \$name | cut -d '_' -f 1)
        taxid=\$(echo \$name | cut -d '_' -f 2)

        sed "s/^>.*/>\${taxid}/" \$sequence >> sco_\${sco}.faa
    done
    """

    stub:
    """
    touch sco_0000.faa
    """
}