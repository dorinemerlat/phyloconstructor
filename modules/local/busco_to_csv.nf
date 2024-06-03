process BUSCO_TO_CSV {
    tag "${id}_${dataset}_${mode}"
    label 'jq'
    scratch false

    input:
    tuple val(id), val(dataset), val(mode), path(json)

    output:
    tuple val(dataset), val(mode), path("busco_${id}_${dataset}_${mode}.csv")

    script:
    """
    if grep -q '"Complete percentage":' $json; then
        busco_res=\$(jq -r '.results | [."Complete percentage", ."Single copy percentage", ."Multi copy percentage", ."Fragmented percentage", ."Missing percentage", .n_markers] | @csv' $json )
    else
        busco_res=\$(jq -r '.results | [.Complete, ."Single copy", ."Multi copy", .Fragmented, .Missing, .n_markers] | @csv' $json )
    fi 
    
    echo \$busco_res | awk -v id=$id -v mode=$mode -v dataset=$dataset 'BEGIN {print "id,mode,dataset,complete,single_copy,multi_copyfragmented,missing,n_markers"} {print id "," mode "," dataset "," \$0}' \
            > busco_${id}_${dataset}_${mode}.csv
    """

    stub:
    """
    touch busco_${id}_${dataset}_${mode}.csv
    """
}