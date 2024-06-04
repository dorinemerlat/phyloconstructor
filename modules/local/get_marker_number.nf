process GET_MARKER_NUMBER {
    tag "${id}"

    input:
    tuple val(dataset), path(busco_dataset) 

    output:
    tuple val(id), path("marker.csv")

    script:
    """
    wc -l ${busco_dataset}/busco_downloads/lineages/${dataset}/links_to_ODB10.txt > marker.csv
    """

    stub:
    """
    touch marker.csv
    """
}