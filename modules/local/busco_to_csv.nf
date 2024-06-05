process BUSCO_TO_CSV {
    tag "${type}"
    label 'phyloconstructor'
    scratch false

    input:
    tuple val(taxid), val(id), path(json), val(type)

    output:
    tuple val(id), path("busco_${id}.csv")

    script:
    """
    if grep -q '"Complete percentage":' $json; then
        busco_res=\$(jq -r '.results | [."Complete percentage", ."Single copy percentage", ."Multi copy percentage", ."Fragmented percentage", ."Missing percentage", .n_markers] | @csv' $json )
    else
        busco_res=\$(jq -r '.results | [.Complete, ."Single copy", ."Multi copy", .Fragmented, .Missing, .n_markers] | @csv' $json )
    fi 
    
    echo \$busco_res | awk -v id=$id -v type=$type 'BEGIN {print "id,type_data,complete,single_copy,multi_copyfragmented,missing,n_markers"} {print id "," type "," \$0}' \
            > busco_${id}.csv
    """

    stub:
    """
    touch busco_${id}.csv
    """
}