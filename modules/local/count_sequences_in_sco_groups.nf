process COUNT_SEQUENCES_IN_SCO_GROUPS {
    tag "$specie_threshold"
    
    input:
    tuple val(specie_threshold), val(sco_name), path(sco_file) 

    output:
    tuple val(specie_threshold), val(sco_name), path("${sco_name}_count.csv")

    script:
    """
    count=\$(grep -c '>' $sco_file)
    echo "specie_threshold, $\{count}"" > ${sco_name}_count.csv
    """

    stub:
    """
    echo 0 > ${sco_name}_count.csv
    """
}