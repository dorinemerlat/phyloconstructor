process GATHER_BUSCO_CSV {
    input:
    path(inputs)

    output:
    path("all_busco.csv")

    script:
    """
    first_input=\$(echo "$inputs" | cut -f 1 -d ' ')
    head -n 1 \$first_input > all_busco.csv
    tail -q -n +2 $inputs >> all_busco.csv
    """

    stub:
    """
    touch ${all_busco.csv}
    """
}