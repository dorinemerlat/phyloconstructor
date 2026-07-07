class Utils {

    /*
     * Parse TSV-like ID files produced by FETCH_* processes.
     *
     * Expected columns after the header:
     *   taxid    species_name    data_id
     *
     * The parser:
     *   - skips the header
     *   - ignores empty lines
     *   - skips malformed records
     *   - normalizes species names for downstream use
     *   - removes duplicates
     *   - sorts records to ensure deterministic execution
     *
     * Output tuple format:
     *   taxid, species, data_id
     */
    static def parse_ids_file(ch) {
        ch
            .flatMap { taxid, ids_file ->
                ids_file.text
                    .readLines()
                    .drop(1)
                    .findAll { it.trim() }
                    .collect { line ->

                        def cols = line.trim().split(/\s+/)

                        if (cols.size() < 3) {
                            log.warn "Skipping malformed line in ${ids_file}: ${line}"
                            return null
                        }

                        def organism_taxid = cols[0].trim()

                        def organism_name = cols[1]
                            .trim()
                            .replaceAll(/\./, '-')
                            .replaceAll(/\s+/, '-')
                            .toLowerCase()

                        def data_id = cols[2].trim()

                        [ organism_taxid, organism_name, data_id ]
                    }
                    .findAll { it != null }
            }
            .distinct()
            .toSortedList { a, b ->
                a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2]
            }
            .flatMap { it }
    }

    /*
     * Prefer RefSeq assemblies when available.
     *
     * For each species:
     *   - keep all GCF accessions if at least one exists
     *   - otherwise keep all GCA accessions
     *
     * This mirrors the common practice of prioritizing RefSeq
     * assemblies over GenBank assemblies whenever possible.
     *
     * Input tuple format:
     *   taxid, species, accession
     *
     * Output tuple format:
     *   taxid, species, accession
     */
    static def keep_refseq_if_available(ch) {
        ch
            .groupTuple(by: [0, 1])
            .map { taxid, organism_name, data_ids ->

                def ids = data_ids
                    .collect { it.toString().trim() }
                    .unique()
                    .sort()

                def gcf_ids = ids.findAll { it.startsWith("GCF_") }

                def selected_ids = gcf_ids ?
                    gcf_ids :
                    ids.findAll { it.startsWith("GCA_") }

                selected_ids.collect { data_id ->
                    [ taxid, organism_name, data_id ]
                }
            }
            .flatMap { it }
    }

    /*
     * Add a source label to transcriptome-derived datasets.
     *
     * Input tuple format:
     *   taxid, species, data_id, path
     *
     * Output tuple format:
     *   taxid, species, data_id, path, source
     */
    static def add_source(ch, source) {
        ch.map { taxid, specie, id, path ->
            [ taxid, specie, id, path, source ]
        }
    }

    /*
     * Add metadata required by BUSCO.
     *
     * Input tuple format:
     *   taxid, species, data_id, path
     *
     * Output tuple format:
     *   taxid, species, data_id, path,
     *   source, busco_mode, busco_dataset
     */
    static def add_source_and_busco_dataset(ch, source, mode, dataset) {
        ch.map { taxid, specie, id, path ->
            [ taxid, specie, id, path, source, mode, dataset ]
        }
    }

    /*
     * Group BUSCO tables by species.
     *
     * BUSCO results originating from multiple data sources
     * (e.g. UniProt, NCBI, transcriptomes) are grouped together
     * for downstream selection of the best dataset.
     *
     * Tables are sorted by file name to ensure deterministic
     * behaviour across runs.
     *
     * Output tuple format:
     *   label, taxid, species, [tables]
     */
    static def group_busco_tables_by_species(busco_tables, label) {
        busco_tables
            .map { taxid, specie, data_id, table, source ->
                [ taxid, specie, table ]
            }
            .groupTuple(by: 0)
            .map { taxid, species, tables ->

                def specie = species.unique().sort()[0]

                def sorted_tables = tables
                    .unique { it.name }
                    .sort { it.name }

                [ label, taxid, specie, sorted_tables ]
            }
    }
}