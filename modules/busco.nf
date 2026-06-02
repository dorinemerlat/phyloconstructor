process BUSCO {
    tag "${specie}/${source}/${data_id}"
    stageInMode 'copy'

    cpus { mode == 'genome' ? 6 : 10 }
    memory {
        mode == 'genome'
            ? "${100 + (50 * (task.attempt - 1))} GB"
            : "${20 + (10 * (task.attempt - 1))} GB"
    }

    time '10h'

    input:
    tuple val(taxid), val(specie), val(data_id), path(fasta), val(source), val(mode), val(dataset)

    output:
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.json"), val(source), emit: json
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}.txt"), val(source), emit: txt
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}_full_table.tsv"), val(source), emit: full_table
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}_single_copy_busco_sequences.fasta"), val(source), emit: single_copy_busco_sequences
    tuple val(taxid), val(specie), val(data_id), path("${specie}_${source}_${data_id}_multi_copy_busco_sequences.fasta"), val(source), emit: multi_copy_busco_sequences
    
    script:
    def busco_out = "busco_${specie}_${source}_${data_id}_${mode}"
    def single_copy_dir = "${busco_out}/run_${dataset}/busco_sequences/single_copy_busco_sequences"
    def multi_copy_dir = "${busco_out}/run_${dataset}/busco_sequences/multi_copy_busco_sequences"
    def single_copy_file = "${specie}_${source}_${data_id}_single_copy_busco_sequences.fasta"
    def multi_copy_file = "${specie}_${source}_${data_id}_multi_copy_busco_sequences.fasta"
    def prefix = "${specie}.${source}.${data_id}"

    """
    set -euo pipefail

    module load busco

    export TMPDIR="\$PWD/tmp"
    mkdir -p "\$TMPDIR"
    export TMP="\$TMPDIR"
    export TEMP="\$TMPDIR"

    sed 's|/|_|g' ${fasta} > ${specie}_${source}_${data_id}.fasta

    busco \\
        -i ${specie}_${source}_${data_id}.fasta \\
        -o ${busco_out} \\
        -m ${mode} \\
        -l ${dataset} \\
        -c ${task.cpus} \\
        -f \\
        --offline \\
        --download_path /shared/projects/metainvert/phyloconstructor2/data/busco_downloads

    mv ${busco_out}/short_summary*.json ${specie}_${source}_${data_id}.json
    mv ${busco_out}/short_summary*.txt ${specie}_${source}_${data_id}.txt
    mv ${busco_out}/run_${dataset}/full_table.tsv ${specie}_${source}_${data_id}_full_table.tsv

    format_busco_sequences() {
        local busco_dir="\$1"
        local out_file="\$2"
        local mode_value="\$3"

        > "\$out_file"

        if [[ ! -d "\$busco_dir" ]] || ! find "\$busco_dir" -type f -name "*.faa" -print -quit | grep -q .; then
            return 0
        fi

        find "\$busco_dir" -type f -name "*.faa" | sort | while read -r faa; do
            busco_id=\$(basename "\$faa" .faa)

            if [[ "\$mode_value" == "genome" ]]; then
                sed "s|^>|>${prefix}.|g" "\$faa" >> "\$out_file"
            else
                sed "s|^>|>${prefix}.\${busco_id}.|g" "\$faa" >> "\$out_file"
            fi
        done
    }

    format_busco_sequences "${single_copy_dir}" "${single_copy_file}" "${mode}"
    format_busco_sequences "${multi_copy_dir}" "${multi_copy_file}" "${mode}"
    """
}