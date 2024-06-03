process BUSCO {
    tag "${id}_${dataset}_${mode}"
    cpus 30
    memory '50 GB'
    time '2d'
    label 'busco'
    scratch true
    stageOutMode 'move'

    input:
    tuple val(id), path(genome), val(dataset), val(mode), path(busco_downloads)

    output:
    tuple val(id), val(dataset), val(mode), path("busco_${id}_${dataset}_${mode}.json"),     emit: json
    tuple val(id), val(dataset), val(mode), path("busco_${id}_${dataset}_${mode}.txt"),      emit: txt
    tuple val(id), val(dataset), val(mode), path("${id}/run_${dataset}/busco_sequences/"),   emit: busco_sequences
    tuple val(id), val(dataset), val(mode), path("${id}/run_${dataset}/full_table.tsv"),     emit: full_table

    script:
    """
    busco -i $genome -l $dataset -o $id -m $mode -c $task.cpus --offline --download_path $busco_downloads -f

    # copy BUSCO files to output
    cp ${id}/short_summary.*.txt busco_${id}_${dataset}_${mode}.txt
    cp ${id}/short_summary.*.json busco_${id}_${dataset}_${mode}.json

    """

    stub:
    """
    touch busco_${id}_${dataset}_${mode}.json
    touch busco_${id}_${dataset}_${mode}.txt
    """
}