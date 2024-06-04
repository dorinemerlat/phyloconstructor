process BUSCO {
    tag "${id}"
    cpus 30
    memory '50 GB'
    time '2d'
    label 'busco'
    scratch true
    stageOutMode 'move'

    input:
    tuple val(id), path(sequences), val(dataset), val(mode), path(busco_downloads)

    output:
    tuple val(id), path("busco_${id}.json"),     emit: json
    tuple val(id), path("busco_${id}.txt"),      emit: txt
    tuple val(id), path("${id}/run_${dataset}/busco_sequences/"),   emit: busco_sequences

    script:
    """
    busco -i $sequences -l $dataset -o $id -m $mode -c $task.cpus --offline --download_path $busco_downloads -f

    # copy BUSCO files to output
    cp ${id}/short_summary.*.txt busco_${id}.txt
    cp ${id}/short_summary.*.json busco_${id}.json
    """

    stub:
    """
    touch busco_${id}.json
    touch busco_${id}.txt
    """
}