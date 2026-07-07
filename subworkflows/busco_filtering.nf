import Utils

include { BUSCO                          } from '../modules/busco_filtering/busco'
include { FILTER_SINGLE_COPY_BUSCO_HITS  } from '../modules/busco_filtering/filter_single_copy_busco_hits'
include { FILTER_ALL_COMPLETE_BUSCO_HITS } from '../modules/busco_filtering/filter_all_complete_busco_hits'
include { SELECT_BEST_BUSCO_HITS         } from '../modules/busco_filtering/select_best_busco_hits'
include { SELECT_ORTHOGROUPS             } from '../modules/busco_filtering/select_orthogroups'

workflow BUSCO_FILTERING {

    take:
    all_busco_inputs

    main:

    BUSCO(all_busco_inputs)

    FILTER_SINGLE_COPY_BUSCO_HITS(BUSCO.out.full_table)
    FILTER_ALL_COMPLETE_BUSCO_HITS(BUSCO.out.full_table)

    Utils.group_busco_tables_by_species(
        FILTER_ALL_COMPLETE_BUSCO_HITS.out,
        'all_complete'
    )
        .concat(
            Utils.group_busco_tables_by_species(
                FILTER_SINGLE_COPY_BUSCO_HITS.out,
                'only_single_copy'
            )
        )
        .set { complete_busco_tables_by_species }

    SELECT_BEST_BUSCO_HITS(complete_busco_tables_by_species)

    SELECT_BEST_BUSCO_HITS.out
        .map { label, taxid, specie, table ->
            [ label, table ]
        }
        .groupTuple(by: 0)
        .map { label, tables ->
            def sorted_tables = tables.unique { it.name }.sort { it.name }
            [ label, sorted_tables, params.busco_dataset_size ]
        }
        .set { all_busco_results }

    SELECT_ORTHOGROUPS(all_busco_results)

    BUSCO.out.single_copy_busco_sequences
        .map { taxid, specie, data_id, fasta, source ->
            fasta
        }
        .collect()
        .map { files ->
            files.unique { it.name }.sort { it.name }
        }
        .set { single_copy_busco_sequences }

    BUSCO.out.multi_copy_busco_sequences
        .map { taxid, specie, data_id, fasta, source ->
            fasta
        }
        .collect()
        .map { files ->
            files.unique { it.name }.sort { it.name }
        }
        .set { multi_copy_busco_sequences }

    emit:
    selected_orthogroups   = SELECT_ORTHOGROUPS.out.tables
    single_copy_sequences  = single_copy_busco_sequences
    multi_copy_sequences   = multi_copy_busco_sequences
}