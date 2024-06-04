process GET_TAXID {
    tag "${id}"

    input:
    tuple val(id), path(csv) 

    output:
    tuple val(id), path('taxid.txt')

    script:
    """
    grep '${id}' ${csv} | cut -f 1 -d ',' > taxid.txt
    """

    stub:
    """
    echo 0000 > taxid.txt
    """
}