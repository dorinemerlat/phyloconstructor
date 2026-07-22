process RNASPADES {
    tag "${specie}"
    label 'rna_tools'
    cpus 6
    memory { "${50 + (50 * (task.attempt - 1))} GB" }
    time '5d'

    input:
    tuple val(taxid), val(specie), path(reads1), path(reads2)

    output:
    tuple val(taxid), val(specie), path("${specie}_rnaspades_transcripts.fasta")

    script:
    """
    check_and_repair_pair() {
        local r1="\$1"
        local r2="\$2"
        local prefix="\$3"

        local n1_lines
        local n2_lines
        local n1
        local n2

        n1_lines=\$(gzip -cd "\$r1" | wc -l)
        n2_lines=\$(gzip -cd "\$r2" | wc -l)

        if (( n1_lines % 4 != 0 || n2_lines % 4 != 0 ))
        then
            echo "ERROR: FASTQ line count is not divisible by four" >&2
            echo "\$r1: \$n1_lines lines" >&2
            echo "\$r2: \$n2_lines lines" >&2
            exit 1
        fi

        n1=\$((n1_lines / 4))
        n2=\$((n2_lines / 4))

        if [[ "\$n1" -eq "\$n2" ]]
        then
            printf "%s\\n%s\\n" "\$r1" "\$r2"
            return
        fi

        echo "Unbalanced pair detected. Running repair.sh." >&2

        local repaired1="\${prefix}_1.repaired.fastq.gz"
        local repaired2="\${prefix}_2.repaired.fastq.gz"
        local singletons="\${prefix}_singletons.fastq.gz"

        repair.sh \\
            in1="\$r1" \\
            in2="\$r2" \\
            out1="\$repaired1" \\
            out2="\$repaired2" \\
            outs="\$singletons" \\
            repair=t

        if [[ ! -s "\$repaired1" || ! -s "\$repaired2" ]]
        then
            echo "ERROR: repair.sh produced empty paired FASTQ files" >&2
            exit 1
        fi

        printf "%s\\n%s\\n" "\$repaired1" "\$repaired2"
    }

    r1_files=( ${reads1.join(' ')} )
    r2_files=( ${reads2.join(' ')} )

    if [[ "\${#r1_files[@]}" -ne "\${#r2_files[@]}" ]]
    then
        echo "ERROR: different numbers of R1 and R2 files" >&2
        exit 1
    fi

    fixed_r1=()
    fixed_r2=()

    # Validate and repair each sequencing run independently.
    for i in "\${!r1_files[@]}"
    do
        prefix="${specie}_pair\$((i + 1))"

        mapfile -t repaired_pair < <(
            check_and_repair_pair \\
                "\${r1_files[\$i]}" \\
                "\${r2_files[\$i]}" \\
                "\$prefix"
        )

        fixed_r1+=("\${repaired_pair[0]}")
        fixed_r2+=("\${repaired_pair[1]}")
    done

    rnaspades_inputs=()

    for i in "\${!fixed_r1[@]}"
    do
        rnaspades_inputs+=(
            "-1" "\${fixed_r1[\$i]}"
            "-2" "\${fixed_r2[\$i]}"
        )
    done

    # Assemble all paired RNA-seq runs for the species.
    rnaspades.py \\
        "\${rnaspades_inputs[@]}" \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga() as long} \\
        -o "rnaspades_${specie}"

    cp "rnaspades_${specie}/transcripts.fasta" \\
        "${specie}_rnaspades_transcripts.fasta"

    test -s "${specie}_rnaspades_transcripts.fasta" || {
        echo "ERROR: rnaSPAdes produced an empty transcriptome" >&2
        exit 1
    }
    """

    stub:
    """
    command -v rnaspades.py >/dev/null
    command -v repair.sh >/dev/null
    command -v gzip >/dev/null

    touch "${specie}_rnaspades_transcripts.fasta"
    """
}