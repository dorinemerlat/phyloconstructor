```groovy
import Utils
import groovy.json.JsonSlurper

include { BUSCO                          } from '../modules/busco_filtering/busco'
include { FILTER_SINGLE_COPY_BUSCO_HITS  } from '../modules/busco_filtering/filter_single_copy_busco_hits'
include { FILTER_ALL_COMPLETE_BUSCO_HITS } from '../modules/busco_filtering/filter_all_complete_busco_hits'
include { SELECT_BEST_BUSCO_HITS         } from '../modules/busco_filtering/select_best_busco_hits'
include { SELECT_ORTHOGROUPS             } from '../modules/busco_filtering/select_orthogroups'


workflow BUSCO_FILTERING {

    take:
    all_busco_inputs
    busco_filtering_strategy

    main:

    /*
     * Run BUSCO on every candidate dataset.
     */
    BUSCO(all_busco_inputs)


    /*
     * Read the OrthoDB dataset size from the first BUSCO JSON report.
     *
     * All BUSCO analyses use the same lineage dataset, so every JSON
     * report is expected to contain the same number_of_buscos value.
     */
    BUSCO.out.json
        .first()
        .map { taxid, specie, data_id, json_file, source ->

            def json_data = new JsonSlurper().parse(json_file.toFile())

            def find_json_value
            find_json_value = { object, key ->

                if (object instanceof Map) {
                    if (object.containsKey(key)) {
                        return object[key]
                    }

                    for (value in object.values()) {
                        def result = find_json_value(value, key)

                        if (result != null) {
                            return result
                        }
                    }
                }

                if (object instanceof Collection) {
                    for (value in object) {
                        def result = find_json_value(value, key)

                        if (result != null) {
                            return result
                        }
                    }
                }

                return null
            }

            def number_of_buscos = find_json_value(
                json_data,
                'number_of_buscos'
            )

            if (number_of_buscos == null) {
                throw new IllegalStateException(
                    "Could not find 'number_of_buscos' in BUSCO report: ${json_file}"
                )
            }

            def dataset_size

            try {
                dataset_size = number_of_buscos.toString().toInteger()
            } catch (Exception error) {
                throw new IllegalStateException(
                    "Invalid number_of_buscos value '${number_of_buscos}' " +
                    "in BUSCO report: ${json_file}",
                    error
                )
            }

            log.info(
                "BUSCO lineage dataset size: ${dataset_size} orthologs " +
                "(read from ${json_file.name})"
            )

            dataset_size
        }
        .set { busco_dataset_size }


    /*
     * Apply the BUSCO filtering strategy selected by the user.
     *
     * Allowed modes:
     *
     *   only_single_copy
     *       Retain complete single-copy BUSCO hits only.
     *
     *   all_complete
     *       Retain all complete BUSCO hits, including duplicated hits.
     *
     *   both
     *       Run both strategies independently.
     *
     * Parameter validation is performed in main.nf before the workflow
     * starts, so only valid values can reach this subworkflow.
     */
    if (busco_filtering_strategy == 'only_single_copy') {

        FILTER_SINGLE_COPY_BUSCO_HITS(
            BUSCO.out.full_table
        )

        Utils.group_busco_tables_by_species(
            FILTER_SINGLE_COPY_BUSCO_HITS.out,
            'only_single_copy'
        )
            .set { complete_busco_tables_by_species }

    } else if (busco_filtering_strategy == 'all_complete') {

        FILTER_ALL_COMPLETE_BUSCO_HITS(
            BUSCO.out.full_table
        )

        Utils.group_busco_tables_by_species(
            FILTER_ALL_COMPLETE_BUSCO_HITS.out,
            'all_complete'
        )
            .set { complete_busco_tables_by_species }

    } else {

        /*
         * Run both BUSCO filtering strategies and keep them as separate
         * datasets through the downstream analyses.
         */
        FILTER_SINGLE_COPY_BUSCO_HITS(
            BUSCO.out.full_table
        )

        FILTER_ALL_COMPLETE_BUSCO_HITS(
            BUSCO.out.full_table
        )

        Utils.group_busco_tables_by_species(
            FILTER_SINGLE_COPY_BUSCO_HITS.out,
            'only_single_copy'
        )
            .concat(
                Utils.group_busco_tables_by_species(
                    FILTER_ALL_COMPLETE_BUSCO_HITS.out,
                    'all_complete'
                )
            )
            .set { complete_busco_tables_by_species }
    }


    /*
     * Select the best BUSCO dataset for each species independently
     * within each BUSCO filtering strategy.
     */
    SELECT_BEST_BUSCO_HITS(
        complete_busco_tables_by_species
    )


    /*
     * Group the selected BUSCO tables by filtering strategy.
     *
     * Intermediate format:
     *
     *   [label, list_of_tables]
     */
    SELECT_BEST_BUSCO_HITS.out
        .map { label, taxid, specie, table ->
            [label, table]
        }
        .groupTuple(by: 0)
        .map { label, tables ->

            def sorted_tables = tables
                .unique { it.name }
                .sort { it.name }

            [label, sorted_tables]
        }
        .set { grouped_busco_results }


    /*
     * Add the lineage dataset size extracted from the BUSCO JSON report.
     *
     * Final format:
     *
     *   [label, list_of_tables, BUSCO_dataset_size]
     */
    grouped_busco_results
        .combine(busco_dataset_size)
        .map { label, tables, dataset_size ->
            [label, tables, dataset_size]
        }
        .set { all_busco_results }


    /*
     * Select orthogroups according to the configured BUSCO completeness
     * and gene occupancy thresholds.
     */
    SELECT_ORTHOGROUPS(
        all_busco_results
    )


    /*
     * Collect complete single-copy BUSCO sequences from all candidate
     * datasets. These files are required later to extract sequences for
     * orthogroups selected with the single-copy strategy.
     */
    BUSCO.out.single_copy_busco_sequences
        .map { taxid, specie, data_id, fasta, source ->
            fasta
        }
        .collect()
        .map { files ->
            files
                .unique { it.name }
                .sort { it.name }
        }
        .set { single_copy_busco_sequences }


    /*
     * Collect duplicated BUSCO sequences from all candidate datasets.
     * These files are required when the all_complete strategy is used.
     */
    BUSCO.out.multi_copy_busco_sequences
        .map { taxid, specie, data_id, fasta, source ->
            fasta
        }
        .collect()
        .map { files ->
            files
                .unique { it.name }
                .sort { it.name }
        }
        .set { multi_copy_busco_sequences }


    emit:
    selected_orthogroups  = SELECT_ORTHOGROUPS.out.tables
    single_copy_sequences = single_copy_busco_sequences
    multi_copy_sequences  = multi_copy_busco_sequences
}
```
