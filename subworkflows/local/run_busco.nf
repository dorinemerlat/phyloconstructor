/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUSCO             } from "../../modules/local/busco"
include { BUSCO_TO_CSV      } from "../../modules/local/busco_to_csv"
include { GATHER_BUSCO_CSV  } from "../../modules/local/gather_busco_csv"


def join_taxid_mode(taxid, busco) {
    return taxid.join(busco).map {id , taxId, sequences, mode  -> [taxId, id, sequences, mode]}
}


workflow RUN_BUSCO {
    take:
        busco_channel
        busco_datasets
        marker_number
        taxid

    main:
        // run busco
        BUSCO(busco_channel)
        
        // set the taxid as id for each item
        join_taxid_mode(taxid, BUSCO.out.json)
            .set { all_busco_json }
        
        // concert the result to csv and gather them
        BUSCO_TO_CSV(all_busco_json)
        GATHER_BUSCO_CSV(BUSCO_TO_CSV.out.map {id, csv -> csv}.toList())

        join_taxid_mode(taxid, BUSCO.out.busco_sequences)
            .map { taxid, id, busco, mode -> [taxid, busco] }
            .groupTuple()
            .combine(marker_number)
            .set { all_busco_sequences }

    // emit:
    //     busco_channel = busco_channel
    //     busco_sequences = all_busco_sequences
}
