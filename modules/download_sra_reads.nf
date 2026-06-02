process DOWNLOAD_SRA_READS {
    tag "$specie"
    cache 'lenient'

    memory { "${20 + (50 * (task.attempt - 1))} GB" }
    cpus 4
    maxForks 20

    input:
    tuple val(taxid), val(specie), val(sra_id)

    output:
    tuple val(taxid), val(specie), val(sra_id), path("${specie}_${sra_id}_1.fastq.gz"), path("${specie}_${sra_id}_2.fastq.gz"), emit: "sra_reads"
    tuple val(taxid), val(specie), val(sra_id), path("${specie}_${sra_id}.read_counts.tsv"), emit: "read_counts"

    script:
    """
    module load sra-tools/3.4.1 pigz bbmap

    count_reads() {
        local fq="\$1"

        if [[ ! -f "\$fq" ]]; then
            echo 0
            return
        fi

        local n_lines
        n_lines=\$(zcat "\$fq" | wc -l)

       if (( n_lines % 4 != 0 )); then
            echo "ERROR: FASTQ line count is not divisible by 4: \$fq" >&2
            echo "\$fq: \$n_lines lines" >&2
            exit 1
        fi

        echo \$((n_lines / 4))
    }

    # export TMPDIR=\$PWD/tmp_sra_${sra_id}
    # export TEMP="\$TMPDIR"
    # export TMP="\$TMPDIR"

    # mkdir -p "\$TMPDIR"

    # prefetch -f ALL ${sra_id} --max-size 200G --output-directory .

    # fasterq-dump \\
    #     --split-3 \\
    #     --threads ${task.cpus} \\
    #     --mem 1000M \\
    #     --disk-limit-tmp 500G \\
    #     --disk-limit 500G \\
    #     -x \\
    #     --temp "\$TMPDIR" \\
    #     ${sra_id} \\
    #     -O . \\
    #     -o ${specie}_${sra_id}

    # pigz -p ${task.cpus} *.fastq

    cd /shared/projects/metainvert/phyloconstructor2/cache/download_sra_reads/${specie}

    r1="${specie}_${sra_id}_1.fastq.gz"
    r2="${specie}_${sra_id}_2.fastq.gz"

    if [[ ! -f "\$r1" || ! -f "\$r2" ]]; then
        echo "ERROR: expected paired-end files were not produced" >&2
        echo "Missing or absent: \$r1 / \$r2" >&2
        ls -lh >&2
        exit 1
    fi

    n1=\$(count_reads "\$r1")
    n2=\$(count_reads "\$r2")

    echo "Initial read counts:"
    echo "  \$r1: \$n1"
    echo "  \$r2: \$n2"

    if [[ "\$n1" -ne "\$n2" ]]; then
        echo "Unbalanced paired-end reads detected. Running repair.sh..."

        repaired1="${specie}_${sra_id}_1.repaired.fastq.gz"
        repaired2="${specie}_${sra_id}_2.repaired.fastq.gz"

        repair.sh \\
            in1="\$r1" \\
            in2="\$r2" \\
            out1="\$repaired1" \\
            out2="\$repaired2" \\
            repair=t \\
            ain=t

        if [[ ! -s "\$repaired1" || ! -s "\$repaired2" ]]; then
            echo "ERROR: repair.sh produced empty repaired files" >&2
            echo "This run may not be true paired-end." >&2
            echo "R1 reads before repair: \$n1" >&2
            echo "R2 reads before repair: \$n2" >&2
            echo -e "${taxid}\\t${specie}\\t${sra_id}\\t\$n1\\t\$n2" > "${specie}_${sra_id}.read_counts.tsv"
        fi

        # rm -f "\$r1" "\$r2"
        mv "\$r1" "${specie}_${sra_id}_1.unrepaired.fastq.gz"
        mv "\$r2" "${specie}_${sra_id}_2.unrepaired.fastq.gz"

        mv "\$repaired1" "\$r1"
        mv "\$repaired2" "\$r2"

        n1=\$(count_reads "\$r1")
        n2=\$(count_reads "\$r2")

        echo "Read counts after repair:"
        echo "  \$r1: \$n1"
        echo "  \$r2: \$n2"

        if [[ "\$n1" -ne "\$n2" ]]; then
            echo "ERROR: reads are still unbalanced after repair" >&2
            echo -e "${taxid}\\t${specie}\\t${sra_id}\\t\$n1\\t\$n2" > "${specie}_${sra_id}.read_counts.tsv"
        fi
    else
        echo "Paired-end reads are balanced. No repair needed."
    fi

    echo -e "${taxid}\\t${specie}\\t${sra_id}\\t\$n1\\t\$n2" > "${specie}_${sra_id}.read_counts.tsv"

    # rm -rf ${sra_id}
    # rm -rf "\$TMPDIR"
    """
}