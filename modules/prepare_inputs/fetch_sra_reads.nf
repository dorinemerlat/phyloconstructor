process FETCH_SRA_READS {
    tag "${taxid}"
    label 'entrez_direct'
    memory '2 GB'
    time '2h'

    input:
    val taxid

    output:
    tuple val(taxid), path("${taxid}.sra_runs.tsv"), emit: ids
    tuple val(taxid), path("${taxid}.sra_runs.out"), emit: info

    script:
    """
    search_sra() {
        local query_taxid="\$1"

        esearch -db sra \\
            -query '((((txid'"\${query_taxid}"'[Organism:exp]) AND "paired"[Layout]) AND "illumina"[Platform]) AND "rna data"[Filter]) AND "filetype fastq"[Properties]' \\
            > esearch.out
    }

    # Retry the Entrez search once after a randomized delay.
    search_sra "${taxid}" || {
        sleep \$(shuf -i 5-30 -n 1)
        search_sra "${taxid}"
    }

    count=\$(grep -o '<Count>[0-9]*</Count>' esearch.out | head -n 1 | sed 's/<[^>]*>//g')
    count=\${count:-0}

    if [[ "\$count" != "0" ]]
    then
        efetch -format runinfo < esearch.out \\
            > "${taxid}.sra_runs.out"

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
        ' "${taxid}.sra_runs.out" \\
            > "${taxid}.sra_runs.tsv"
    else
        printf "taxid\\tspecie\\tsra\\n" \\
            > "${taxid}.sra_runs.tsv"

        touch "${taxid}.sra_runs.out"
    fi
    """

    stub:
    """
    command -v esearch >/dev/null
    command -v efetch >/dev/null

    touch "${taxid}.sra_runs.tsv"
    touch "${taxid}.sra_runs.out"
    """
}