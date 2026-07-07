process FETCH_SRA_READS {
    tag ""

    input:
    val(taxid)

    output:
    tuple val(taxid), path("${taxid}.sra_runs.tsv"), emit: ids
    tuple val(taxid), path("${taxid}.sra_runs.out"), emit: info

    script:
    """
    module load entrez-direct/22.4

    search.sra() {
        local taxid="\$1"

        esearch -db sra \\
            -query '((((txid'"\\\${taxid}"'[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \\
            > esearch.out
    }

    {
        search.sra ${taxid}
    } || {
        sleep \$(shuf -i 5-30 -n 1)
        search.sra ${taxid}
    }

    count=\$(grep "<Count>" esearch.out | cut -d '>' -f 2 | cut -d '<' -f 1)

    if [[ "\$count" != "0" ]]; then

        efetch -format runinfo < esearch.out > ${taxid}.sra_runs.out

        awk -F ',' '
            BEGIN { OFS="\\t" }

            NR == 1 {
                print "taxid","specie","sra"
                next
            }

            {
                run_accession    = \$1
                spots_with_mates = \$6
                layout           = \$16
                taxid            = \$28
                specie           = \$29

                if (layout != "PAIRED")
                    next

                if (spots_with_mates == "" || spots_with_mates == "0")
                    next

                specie = tolower(specie)
                gsub(/ /, "-", specie)
                gsub(/\\./, "-", specie)

                print taxid, specie, run_accession
            }

        ' ${taxid}.sra_runs.out > ${taxid}.sra_runs.tsv

    else

        echo -e "taxid\\tspecie\\tsra" > ${taxid}.sra_runs.tsv

    fi
    """
}