/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GROUPING_ORTHOLOGS_BY_SPECIES             } from "../../modules/local/grouping_orthologs_by_species"
include { GROUPING_ORTHOLOGS_BY_SCO_GROUPS          } from "../../modules/local/grouping_orthologs_by_sco_groups"
include { COUNT_SEQUENCES_IN_SCO_GROUPS             } from "../../modules/local/count_sequences_in_sco_groups"


workflow FILTER_SCO_SEQUENCES {
    take:
        busco_sequences
        threshold_species
        threshold_orthogroups

    main:
        GROUPING_ORTHOLOGS_BY_SPECIES(busco_sequences)
        
        threshold_orthogroups = Channel.from(0..10).map { it*10 }
        threshold_species = threshold_orthogroups

        // selection of species with complete percentage of single copy orthologs >= threshold
        GROUPING_ORTHOLOGS_BY_SPECIES.out
            .map { taxid, proteins, csv -> [ taxid, proteins, csv.readLines()[1].split(',') ] }
            .map { taxid, proteins, csv -> [ taxid, proteins, csv[1].toInteger(), csv[0].toInteger() ] } // csv[0] = sco_count, csv[1] = sco_percentage
            .filter { it[3] != 0 } // remove specie where sco == 0 
            .combine(threshold_species)
            .filter { it[2] >= it[4] } // filter species to keep (complete percentage of single copy orthologs >= threshold)
            .map { taxid, proteins, sco_perc, sco_count, threshold_specie -> [threshold_specie, proteins] }
            .groupTuple()
            .set { sco_for_grouping_by_specie }

        // concatenate all the single copy orthologs for the selected species
        GROUPING_ORTHOLOGS_BY_SCO_GROUPS(sco_for_grouping_by_specie)

        // filter to keep only those with a percentage of sequences in each ortholog groups greater than a given threshold
        GROUPING_ORTHOLOGS_BY_SCO_GROUPS.out
            .transpose()
            .map { [threshold_specie, sco_group_dir] -> [threshold_specie, sco_group_dir.getBaseName(), sco_group_dir] }
            .set { sco_groups }

        COUNT_SEQUENCES_IN_SCO_GROUPS(sco_groups)

        sco_groups
            .join(COUNT_SEQUENCES_IN_SCO_GROUPS.out, by: [0,1])
            .map { threshold_specie, sco_name, sco_dir, sco_count -> [threshold_specie, sco_name, sco_dir, sco_count.readLines()[0].split(",")[1].toInteger()] }
            .set { sco_groups_with_count }

        sco_groups_with_count
            .map { threshold_specie, sco_name, sco_dir, sco_count -> threshold_specie, sco_count }
            .max { it[1] }
            .set { sco_max_count }

        sco_groups_with_count
            .join(sco_max_count)
            .map { threshold_specie, sco_name, sco_dir, sco_count, sco_max -> [threshold_specie, sco_name, sco_dir, sco_count*100/sco_max] }
            .set { sco_groups_with_percentage }

        // save the sco_count for each threshold_specie
        sco_groups_with_percentage
            .filter { it[2] >= threshold_orthogroups }
            .map { name, dir, percentage -> ["threshold_specie_${threshold_species}_threshold_orthogroup_${threshold_orthogroups}", name, dir] }
            .set { good_sco_groups }

        good_sco_groups.view()

    emit:
        good_sco_groups = good_sco_groups
}
