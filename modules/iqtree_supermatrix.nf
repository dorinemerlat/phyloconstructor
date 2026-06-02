process IQTREE_SUPERMATRIX {
    tag "${job_name}"
    cpus 20
    memory { "${50 * task.attempt} GB" }
    scratch false
    time '1d'

    input:
    tuple val(job_name), path(aln), path(partition), val(model_outgroups)

    output:
    tuple val(job_name), path("${outgroup}.treefile") 
    tuple val(job_name), path("${outgroup}.*") 

    script:
    """
    module load iqtree

    awk '{print "LG, " \$0}' ${partition} > partition_iqtree.txt

    present_outgroups=\$(echo "${model_outgroups}" | tr ',' '\\n' | while read og; do
        if grep -q "^>\${og}\$\\|^>\${og}[[:space:]]" ${aln}; then
            echo "\${og}"
        fi
    done | paste -sd "," -)

    if [[ -n "\${present_outgroups}" ]]; then
        outgroup_option="-o \${present_outgroups}"
    else
        outgroup_option=""
    fi

    echo "Outgroups requested: ${model_outgroups}"
    echo "Outgroups present: \${present_outgroups}"

    iqtree3 \\
        -s ${aln} \\
        -p partition_iqtree.txt \\
        -m MFP+MERGE \\
        -B 1000 \\
        --alrt 1000 \\
        --bnni \\
        --prefix ${job_name}_supermatrix \\
        -T ${task.cpus} \\
        \${outgroup_option}
    """
}