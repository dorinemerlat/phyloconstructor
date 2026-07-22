process PHYKIT_CONCAT {
    tag "${label}/${job_name}"
    label 'sequence_tools'
    memory { "${10 + (5 * (task.attempt - 1))} GB" }
    time '6h'

    input:
    tuple val(label), val(job_name), path(aln, stageAs: "input/*")

    output:
    tuple val(label), val(job_name), path("${job_name}.fa"), path("${job_name}.partition"), emit: main_output
    tuple val(label), val(job_name), path("${job_name}.occupancy"), emit: occupancy

    script:
    """
    mkdir -p cleaned

    # Remove source metadata from sequence identifiers before concatenation.
    for f in input/*
    do
        seqkit replace \\
            -p '^(.*?)_source=.*' \\
            -r '\$1' \\
            "\$f" \\
            > "cleaned/\$(basename "\$f")"
    done

    # Provide PhyKIT with a reproducibly ordered alignment list.
    find cleaned -type f | sort > alignments.list

    phykit create_concat \\
        -a alignments.list \\
        -p "${job_name}"

    # Wrap the concatenated FASTA sequences at 60 characters.
    seqkit seq -w 60 "${job_name}.fa" > "${job_name}.fa.tmp"
    mv "${job_name}.fa.tmp" "${job_name}.fa"

    test -s "${job_name}.fa" || {
        echo "ERROR: ${job_name}.fa is empty" >&2
        exit 1
    }

    test -s "${job_name}.partition" || {
        echo "ERROR: ${job_name}.partition is empty" >&2
        exit 1
    }

    test -s "${job_name}.occupancy" || {
        echo "ERROR: ${job_name}.occupancy is empty" >&2
        exit 1
    }
    """

    stub:
    """
    command -v phykit >/dev/null
    command -v seqkit >/dev/null
    phykit create_concat --help >/dev/null

    touch "${job_name}.fa"
    touch "${job_name}.partition"
    touch "${job_name}.occupancy"
    """
}
