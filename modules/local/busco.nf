process BUSCO {
    tag "${id}"
    cpus 30
    memory '50 GB'
    time '2d'
    label 'busco'
    stageOutMode 'move'

    input:
    tuple val(id), path(sequences), val(dataset), val(mode), path(busco_downloads), val(priority)

    output:
    tuple val(id), path("busco_${id}.json"), val(dataset),                          emit: json
    tuple val(id), path("busco_${id}.txt"),  val(dataset),                          emit: txt
    tuple val(id), path("sco_sequences_priority${priority}_${id}"), val(dataset),   emit: busco_sequences
    script:
    """
    busco -i $sequences -l $dataset -o $id -m $mode -c $task.cpus --offline --download_path $busco_downloads -f

    # copy BUSCO files to output
    cp ${id}/short_summary.*.txt busco_${id}.txt
    cp ${id}/short_summary.*.json busco_${id}.json

    mv ${id}/run_${dataset}/busco_sequences/single_copy_busco_sequences sco_sequences_priority${priority}_${id}
    """

    stub:
    """
    touch busco_${id}.json
    touch busco_${id}.txt
    mkdir single_copy_busco_sequences_${id}
    touch single_copy_busco_sequences_${id}/0000.faa
    """
}