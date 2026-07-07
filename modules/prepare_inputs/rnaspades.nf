process RNASPADES {
    tag "${specie}"
    scratch false
    cpus 6
    stageInMode 'copy'
    memory { "${50 + (50 * (task.attempt - 1))} GB" }
    time '5d'

    input:
    tuple val(taxid), val(specie), path(reads1), path(reads2)

    output:
    tuple val(taxid), val(specie), path("${specie}_rnaspades_transcripts.fasta")

    script:
    """
    module load spades bbmap

    check_and_repair_pair() {
        local r1="\$1"
        local r2="\$2"
        local prefix="\$3"

        echo "Checking pair:"
        echo "  R1: \$r1"
        echo "  R2: \$r2"

        local n1_lines
        local n2_lines
        local n1
        local n2

        n1_lines=\$(zcat "\$r1" | wc -l)
        n2_lines=\$(zcat "\$r2" | wc -l)

        if (( n1_lines % 4 != 0 || n2_lines % 4 != 0 )); then
            echo "ERROR: FASTQ line count is not divisible by 4"
            echo "\$r1: \$n1_lines lines"
            echo "\$r2: \$n2_lines lines"
            exit 1
        fi

        n1=\$((n1_lines / 4))
        n2=\$((n2_lines / 4))

        echo "  R1 reads: \$n1"
        echo "  R2 reads: \$n2"

        if [[ "\$n1" -eq "\$n2" ]]; then
            echo "  Pair is balanced. No repair needed." >&2
            echo "\$r1"
            echo "\$r2"
        else
            echo "  Pair is unbalanced. Running repair.sh..." >&2

            local repaired1="\${prefix}_1.repaired.fastq.gz"
            local repaired2="\${prefix}_2.repaired.fastq.gz"
            local singletons="\${prefix}_singletons.fastq.gz"

            repair.sh \\
                in1="\$r1" \\
                in2="\$r2" \\
                out1="\$repaired1" \\
                out2="\$repaired2" \\
                outs="\$singletons" \\
                repair

            echo "\$repaired1"
            echo "\$repaired2"
        fi
    }

    r1_files=( ${reads1.join(' ')} )
    r2_files=( ${reads2.join(' ')} )

    if [[ "\${#r1_files[@]}" -ne "\${#r2_files[@]}" ]]; then
        echo "ERROR: different number of R1 and R2 files"
        echo "R1 files: \${#r1_files[@]}"
        echo "R2 files: \${#r2_files[@]}"
        exit 1
    fi

    fixed_r1=()
    fixed_r2=()

    for i in "\${!r1_files[@]}"; do
        r1="\${r1_files[\$i]}"
        r2="\${r2_files[\$i]}"

        prefix="${specie}_pair\$((i + 1))"

        mapfile -t repaired_pair < <(
            check_and_repair_pair "\$r1" "\$r2" "\$prefix" | tail -n 2
        )

        fixed_r1+=("\${repaired_pair[0]}")
        fixed_r2+=("\${repaired_pair[1]}")
    done

    rnaspades_inputs=()

    for i in "\${!fixed_r1[@]}"; do
        rnaspades_inputs+=("-1" "\${fixed_r1[\$i]}" "-2" "\${fixed_r2[\$i]}")
    done

    echo "rnaSPAdes inputs:"
    printf '  %s\\n' "\${rnaspades_inputs[@]}"

    rnaspades.py \\
        "\${rnaspades_inputs[@]}" \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        -o rnaspades_${specie}

    cp rnaspades_${specie}/transcripts.fasta \\
       ${specie}_rnaspades_transcripts.fasta
    """
}