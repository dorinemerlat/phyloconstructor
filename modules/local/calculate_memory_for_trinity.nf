process CALCULATE_MEMORY_FOR_TRINITY {
    tag "${id}"

    input:
    tuple val(id), path(fastq1), path(fastq2)

    output:
    tuple val(id), path("memory.txt")

    script:
    """
    # Function to count the number of reads in a FASTQ file
    count_reads() {
        local file=\$1
        local read_count

        # Count the number of lines and divide by 4 to get the number of reads
        read_count=\$(cat "\$file" |wc -l)
        read_count=\$((read_count / 4))
        echo "\$read_count"
    }

    # Count the reads in both FASTQ files
    reads1=\$(count_reads "$fastq1")
    reads2=\$(count_reads "$fastq2")

    # Calculate the total number of reads
    total_reads=\$((reads1 + reads2))
    
    # Calculate the memory needed (1GB per 1M reads)
    memory_needed=\$((total_reads / 1000000))
    memory_needed=\$((memory_needed + 10))
    memory_needed=\$(printf '%.0f\\n' \$memory_needed)

    # Ensure at least 1GB of memory is allocated
    if [ "\$memory_needed" -eq 0 ]; then
        memory_needed=1
    fi

    # Write the memory needed to a file
    echo "\$memory_needed" > memory.txt
    """

    stub:
    """
    echo "1" > memory.txt
    """
}