include { PREPARE_INPUTS  } from '../subworkflows/prepare_inputs'
include { BUSCO_FILTERING } from '../subworkflows/busco_filtering'
include { PHYLOGENY       } from '../subworkflows/phylogeny'


workflow PHYLOCONSTRUCTOR {
    
    take:
    busco_filtering_strategy

    main:

    /*
    * Collect user-provided data and optionally download public
    * proteomes, assemblies, transcriptomes and sequencing reads.
    */
    PREPARE_INPUTS()


    /*
    * Run BUSCO, select the best dataset for each species and build
    * orthogroup tables for the configured completeness and occupancy
    * thresholds.
    */
    BUSCO_FILTERING(
        PREPARE_INPUTS.out.all_busco_inputs,
        busco_filtering_strategy
    )


    /*
    * Infer concatenation-based and coalescent-based phylogenies from
    * the selected BUSCO orthogroups.
    */
    PHYLOGENY(
        BUSCO_FILTERING.out.selected_orthogroups,
        BUSCO_FILTERING.out.single_copy_sequences,
        BUSCO_FILTERING.out.multi_copy_sequences
    )

}