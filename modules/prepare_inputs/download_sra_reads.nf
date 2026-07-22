process DOWNLOAD_SRA_READS {
    tag "${specie}/${sra_id}"
    label 'rna_tools'
    cpus 4
    memory { "${20 + (50 * (task.attempt - 1))} GB" }
    time '1d'
    maxForks 20

    input:
    tuple val(taxid), val(specie), val(sra_id)

    output:
    tuple val(taxid), val(specie), val(sra_id), path("${specie}_${sra_id}_1.fastq.gz"), path("${specie}_${sra_id}_2.fastq.gz"), emit: sra_reads
    tuple val(taxid), val(specie), val(sra_id), path("${specie}_${sra_id}.read_counts.tsv"), emit: read_counts

    script:
    """
    count_reads() {
        local fastq="\$1"

        if [[ ! -f "\$fastq" ]]
        then
            echo 0
            return
        fi

        local line_count
        line_count=\$(gzip -cd "\$fastq" | wc -l)

        if (( line_count % 4 != 0 ))
        then
            echo "ERROR: FASTQ line count is not divisible by four: \$fastq" >&2
            echo "\$fastq: \$line_count lines" >&2
            exit 1
        fi

        echo \$((line_count / 4))
    }

    export TMPDIR="\$PWD/tmp_sra_${sra_id}"
    export TEMP="\$TMPDIR"
    export TMP="\$TMPDIR"

    mkdir -p "\$TMPDIR"

    # Download the SRA archive before converting it to FASTQ.
    prefetch -f ALL "${sra_id}" \\
        --max-size 200G \\
        --output-directory .

    fasterq-dump \\
        --split-3 \\
        --threads ${task.cpus} \\
        --mem 1000M \\
        --disk-limit-tmp 500G \\
        --disk-limit 500G \\
        -x \\
        --temp "\$TMPDIR" \\
        "${sra_id}" \\
        -O . \\
        -o "${specie}_${sra_id}"

    # Compress the generated FASTQ files concurrently.
    for fastq in *.fastq
    do
        gzip -f "\$fastq" &
    done
    wait

    r1="${specie}_${sra_id}_1.fastq.gz"
    r2="${specie}_${sra_id}_2.fastq.gz"

    if [[ ! -f "\$r1" || ! -f "\$r2" ]]
    then
        echo "ERROR: expected paired-end files were not produced" >&2
        echo "Expected: \$r1 and \$r2" >&2
        ls -lh >&2
        exit 1
    fi

    n1=\$(count_reads "\$r1")
    n2=\$(count_reads "\$r2")

    if [[ "\$n1" -ne "\$n2" ]]
    then
        echo "Unbalanced paired-end reads detected. Running repair.sh..."

        unrepaired1="${specie}_${sra_id}_1.unrepaired.fastq.gz"
        unrepaired2="${specie}_${sra_id}_2.unrepaired.fastq.gz"
        repaired1="${specie}_${sra_id}_1.repaired.fastq.gz"
        repaired2="${specie}_${sra_id}_2.repaired.fastq.gz"
        singletons="${specie}_${sra_id}.singletons.fastq.gz"

        repair.sh \\
            in1="\$r1" \\
            in2="\$r2" \\
            out1="\$repaired1" \\
            out2="\$repaired2" \\
            outs="\$singletons" \\
            repair=t \\
            ain=t

        if [[ ! -s "\$repaired1" || ! -s "\$repaired2" ]]
        then
            echo "ERROR: repair.sh produced empty paired FASTQ files" >&2
            printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \\
                "${taxid}" "${specie}" "${sra_id}" "\$n1" "\$n2" \\
                > "${specie}_${sra_id}.read_counts.tsv"
            exit 1
        fi

        mv "\$r1" "\$unrepaired1"
        mv "\$r2" "\$unrepaired2"
        mv "\$repaired1" "\$r1"
        mv "\$repaired2" "\$r2"

        n1=\$(count_reads "\$r1")
        n2=\$(count_reads "\$r2")

        if [[ "\$n1" -ne "\$n2" ]]
        then
            echo "ERROR: reads remain unbalanced after repair" >&2
            printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \\
                "${taxid}" "${specie}" "${sra_id}" "\$n1" "\$n2" \\
                > "${specie}_${sra_id}.read_counts.tsv"
            exit 1
        fi
    fi

    printf "%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "${taxid}" "${specie}" "${sra_id}" "\$n1" "\$n2" \\
        > "${specie}_${sra_id}.read_counts.tsv"

    rm -rf "${sra_id}" "\$TMPDIR"
    """

    stub:
    """
    command -v prefetch >/dev/null
    command -v fasterq-dump >/dev/null
    command -v gzip >/dev/null
    command -v repair.sh >/dev/null

    touch "${specie}_${sra_id}_1.fastq.gz"
    touch "${specie}_${sra_id}_2.fastq.gz"
    touch "${specie}_${sra_id}.read_counts.tsv"
    """
}