include { PREPARE_INPUTS  } from '../subworkflows/prepare_inputs'
include { BUSCO_FILTERING } from '../subworkflows/busco_filtering'
include { PHYLOGENY       } from '../subworkflows/phylogeny'

workflow PHYLOCONSTRUCTOR {

    main:

    PREPARE_INPUTS()

    BUSCO_FILTERING(PREPARE_INPUTS.out.all_busco_inputs)

    PHYLOGENY(
        BUSCO_FILTERING.out.selected_orthogroups,
        BUSCO_FILTERING.out.single_copy_sequences,
        BUSCO_FILTERING.out.multi_copy_sequences
    )

    // emit:
    // busco_results = BUSCO_FILTERING.out.busco_results
}