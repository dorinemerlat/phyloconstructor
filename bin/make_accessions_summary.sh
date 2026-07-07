#!/usr/bin/env bash
set -euo pipefail

output="${1:-data_accessions_summary.tsv}"
tmpdir="summary_tmp"

rm -rf "$tmpdir"
mkdir -p "$tmpdir"

echo -e "taxid\tspecie\tresource_type\tsource\taccession" > "$output"

safe_find() {
    local pattern="$1"
    shift

    for dir in "$@"; do
        [ -d "$dir" ] || continue
        find "$dir" -type f -name "$pattern"
    done
}

count_lines() {
    local file="$1"
    [ -s "$file" ] && wc -l < "$file" || echo 0
}

count_species() {
    local file="$1"
    [ -s "$file" ] && cut -f2 "$file" | sort -u | wc -l || echo 0
}

count_taxids() {
    local file="$1"
    [ -s "$file" ] && awk -F'\t' '$1 != "NA" {print $1}' "$file" | sort -u | wc -l || echo 0
}

report_count() {
    local label="$1"
    local file="$2"

    echo "Found $(count_lines "$file") accessions in $label"
    echo "  species: $(count_species "$file")"
}

# -------------------------
# NCBI genomes
# -------------------------
safe_find "*.ncbi_accessions.tsv" \
    cache/prepare_inputs/fetch_ncbi_assemblies \
    cache/fetch_ncbi_assemblies \
    | while read -r f; do
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 3 {
                taxid = $1
                specie = $2
                accession = $3

                line[taxid, accession] = taxid OFS specie OFS "genome" OFS "NCBI" OFS accession

                if (accession ~ /^GCF_/) {
                    has_gcf[taxid] = 1
                }
            }
            END {
                for (k in line) {
                    split(k, parts, SUBSEP)
                    taxid = parts[1]
                    accession = parts[2]

                    if (has_gcf[taxid]) {
                        if (accession ~ /^GCF_/) print line[k]
                    } else {
                        print line[k]
                    }
                }
            }
        ' "$f"
    done \
    | sort -u > "$tmpdir/ncbi_genomes.tsv"

report_count "NCBI genomes" "$tmpdir/ncbi_genomes.tsv"

# -------------------------
# NCBI proteomes
# -------------------------
safe_find "*.ncbi_accessions_for_proteins.tsv" \
    cache/prepare_inputs/fetch_ncbi_assemblies \
    cache/fetch_ncbi_assemblies \
    | while read -r f; do
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 3 {
                print $1, $2, "proteome", "NCBI", $3
            }
        ' "$f"
    done \
    | sort -u > "$tmpdir/ncbi_proteomes.tsv"

report_count "NCBI proteomes" "$tmpdir/ncbi_proteomes.tsv"

# -------------------------
# UniProt proteomes
# -------------------------
safe_find "*.upids.tsv" \
    cache/prepare_inputs/fetch_uniprot_proteomes \
    cache/fetch_uniprot_proteomes \
    | while read -r f; do
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 3 {
                print $1, $2, "proteome", "UniProt", $3
            }
        ' "$f"
    done \
    | sort -u > "$tmpdir/uniprot.tsv"

report_count "UniProt proteomes" "$tmpdir/uniprot.tsv"

# -------------------------
# TSA transcriptomes
# -------------------------
safe_find "*.tsa_ids.tsv" \
    cache/prepare_inputs/fetch_tsa_transcriptomes \
    cache/fetch_tsa_transcriptomes \
    | while read -r f; do
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 3 {
                print $1, $2, "transcriptome", "TSA", $3
            }
        ' "$f"
    done \
    | sort -u > "$tmpdir/tsa.tsv"

report_count "TSA transcriptomes" "$tmpdir/tsa.tsv"

# -------------------------
# SRA reads
# -------------------------
safe_find "*.sra_runs.tsv" \
    cache/prepare_inputs/fetch_sra_reads \
    cache/fetch_sra_reads \
    | while read -r f; do
        awk -F'\t' '
            BEGIN { OFS="\t" }
            NR == 1 { next }
            NF >= 3 {
                print $1, $2, "reads", "SRA", $3
            }
        ' "$f"
    done \
    | sort -u > "$tmpdir/sra.tsv"

report_count "SRA reads" "$tmpdir/sra.tsv"

# -------------------------
# SRA assembled transcriptomes
# -------------------------
find cache/prepare_inputs/rnaspades \
    -mindepth 1 -maxdepth 1 -type d 2>/dev/null \
    | awk -F'/' '
        BEGIN { OFS="\t" }
        {
            specie = $NF

            if ( specie == "logs" || specie == "log_rnaspades" || specie ~ /failed/ ) {
                next
            }

            print "NA", specie, "transcriptome", "RNAspades", "assembled"
        }
    ' \
    | sort -u > "$tmpdir/rnaspades.tsv"

report_count "RNAspades assembled transcriptomes" "$tmpdir/rnaspades.tsv"

# -------------------------
# Final output
# -------------------------
cat "$tmpdir"/ncbi_genomes.tsv \
    "$tmpdir"/ncbi_proteomes.tsv \
    "$tmpdir"/uniprot.tsv \
    "$tmpdir"/tsa.tsv \
    "$tmpdir"/sra.tsv \
    | sort -u \
    >> "$output"

echo
echo "Summary:"
echo "  total accessions: $(tail -n +2 "$output" | wc -l)"
echo "  total species:    $(tail -n +2 "$output" | cut -f2 | sort -u | wc -l)"
echo "  total taxids:     $(tail -n +2 "$output" | awk -F'\t' '$1 != "NA" {print $1}' | sort -u | wc -l)"
echo
echo "Wrote: $output"