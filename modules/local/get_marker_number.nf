process GET_MARKER_NUMBER {
    tag "${dataset}"

    input:
    tuple val(dataset), path(busco_dataset) 

    output:
    tuple val(dataset), path("marker.csv")

    script:
    """
    more ${busco_dataset}/lineages/${dataset}/links_to_ODB10.txt  | wc -l > marker.csv
    """

    stub:
    """
    touch 
    """
}