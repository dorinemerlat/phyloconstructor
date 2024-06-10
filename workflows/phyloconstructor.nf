
include { DOWNLOAD_DATA     } from "../subworkflows/local/download_data"
include { RUN_BUSCO         } from "../subworkflows/local/run_busco"
include { MUSCLE            } from "../modules/local/muscle"

params.taxid = "61985"
params.outgroup_uniprot = "7227,6945,6669"
params.outgroup_all = "6850,438506,1519145"
params.sco_by_specie = 60
params.sco_by_group = 70
params.ncbi_api_key = "01679e044efe7ad60b87d93a1bb9085f7b09"
params.busco_dataset = "arthropoda_odb10"

workflow PHYLOCONSTRUCTOR {
    model_outgroups = Channel.from(params.outgroup_uniprot.split(','))
    short_distance_outgroups = Channel.from(params.outgroup_all.split(','))
    taxid = Channel.from(params.taxid)

    all = taxid.concat(model_outgroups, short_distance_outgroups)
    all_without_models = taxid.concat(short_distance_outgroups)
    all_without_models_with_ncbi_key = all_without_models.map { taxid -> [taxid, params.ncbi_api_key] }

    DOWNLOAD_DATA(
        all, 
        all_without_models, 
        all_without_models_with_ncbi_key,
        params.busco_dataset)

    RUN_BUSCO(
        DOWNLOAD_DATA.out.busco_channel,
        DOWNLOAD_DATA.out.busco_datasets,
        DOWNLOAD_DATA.out.marker_number,
        DOWNLOAD_DATA.out.taxid)
    
    // RUN_BUSCO.out.busco_sequences
    //     .combine(Channel.from(0..100).filter { it % 10 == 0 })
    //     .combine(Channel.from(0..100).filter { it % 10 == 0 })
    //     .set { busco_sequences_with_all_thresholds }

    // FILTER_SCO_SEQUENCES(RUN_BUSCO.out.busco_sequences_with_all_thresholds)

    // GROUPING_ORTHOLOGS_BY_SPECIES(all_busco_sequences)
    
    // // selection of species with complete percentage of single copy orthologs >= threshold
    // GROUPING_ORTHOLOGS_BY_SPECIES.out
    //     .map { taxid, proteins, csv -> [ taxid, proteins, csv.readLines()[1].split(',') ] }
    //     .map { taxid, proteins, csv -> [ taxid, proteins, csv[1].toInteger(), csv[0].toInteger() ] }
    //     .filter { it[3] != 0 } // remove specie where sco == 0 
    //     .filter { it[2] >= params.sco_by_specie } // filter species to keep (complete percentage of single copy orthologs >= threshold)
    //     .map { taxid, proteins, percentage, number -> proteins }
    //     .collect()
    //     .map { proteins -> [params.sco_by_specie, proteins] }
    //     .set { good_species }

    // // concatenate all the single copy orthologs for the selected species
    // GROUPING_ORTHOLOGS_BY_SCO_GROUPS(good_species)

    // // filter to keep only those with a percentage of sequences in each ortholog groups greater than a given threshold
    // GROUPING_ORTHOLOGS_BY_SCO_GROUPS.out
    //     .flatten()
    //     .map { sco_group_dir -> [params.sco_by_specie, sco_group_dir.getBaseName(), sco_group_dir] }
    //     .set { sco_groups }

    // COUNT_SEQUENCES_IN_SCO_GROUPS(sco_groups)

    // sco_groups.map{ threshold, name, dir -> [name, dir] } 
    //     .join(COUNT_SEQUENCES_IN_SCO_GROUPS.out)
    //     .map { sco_name, sco_dir, count -> [sco_name, sco_dir, count.readLines()[0].toInteger()] }
    //     .set { sco_groups_with_count }

    // sco_groups_with_count
    //     .map { name, file, count -> count }
    //     .max()
    //     .set { sco_max_count }

    // sco_groups_with_count
    //     .combine(sco_max_count)
    //     .map { name, file, count, max -> [name, file, count*100/max] }
    //     .set { sco_groups_with_percentage }

    // sco_groups_with_percentage
    //     .filter { it[2] >= params.sco_by_group }
    //     .map { name, dir, percentage -> ["threshold_specie_${params.sco_by_specie}_threshold_orthogroup_${params.sco_by_group}", name, dir] }
    //     .set { good_sco_groups }

    // good_sco_groups.view()
    // MUSCLE(good_sco_groups)

    // trimming of the alignments
    // TRIMAL(MUSCLE.out)

    // concatenation of the alignments

    // generation of the matrix

    // tree estimation 
}

