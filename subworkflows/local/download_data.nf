/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CHECK IF FASTA FILES ARE VALIDS, REFORMATE THEM AND CALCULATE THEIR SIZE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_PROTEOMES                        } from "../../modules/local/download_proteomes"
include { DOWNLOAD_REFSEQ_PROTEOMES                 } from "../../modules/local/download_refseq_proteomes"
include { DOWNLOAD_GENOMES                          } from "../../modules/local/download_genomes"
include { DOWNLOAD_TSA                              } from "../../modules/local/download_tsa"
include { DOWNLOAD_SRA                              } from "../../modules/local/download_sra"
include { TRINITY                                   } from "../../modules/local/trinity"
// include { TRIMMOMATIC                               } from "../../modules/local/trimmomatic"
include { CD_HIT_EST as CD_HIT_EST_TSA              } from "../../modules/local/cd_hit_est"
include { CD_HIT_EST as CD_HIT_EST_SRA              } from "../../modules/local/cd_hit_est"
include { REFORMAT_FASTA as REFORMAT_GENOMES        } from "../../modules/local/reformat_fasta"
include { REFORMAT_FASTA as REFORMAT_PROTEOMES      } from "../../modules/local/reformat_fasta"
include { REFORMAT_FASTA as REFORMAT_SRA            } from "../../modules/local/reformat_fasta"
include { REFORMAT_FASTA as REFORMAT_TSA            } from "../../modules/local/reformat_fasta"
include { DOWNLOAD_BUSCO_DATASETS                   } from "../../modules/local/download_busco_datasets"
include { GET_TAXID                                 } from "../../modules/local/get_taxid"
include { GET_MARKER_NUMBER                         } from "../../modules/local/get_marker_number"


def transpose_channel_3_args(channel) {
    return channel
        .map { taxid, csv, fasta -> [taxid, fasta] }
        .transpose()
        .map { taxid, fasta -> [fasta.getSimpleName(), fasta] }
        .filter { it[0] != "0000" }
}


def transpose_channel_4_args(channel) {
    return channel
        .map { taxid, csv, fastq1, fastq2 -> [taxid, fastq1, fastq2] }
        .transpose()
        .map { taxid, fastq1, fastq2 -> [fastq1.getSimpleName().replaceAll("_1", ""), fastq1, fastq2] }
        .filter { it[0] != "0000" }
}


def get_csv_3_args(channel) {
    return channel
        .map { taxid, csv, fasta -> [csv, fasta] }
        .transpose()
        .map { csv, fasta -> [fasta.getSimpleName(), csv] }
        .filter { it[0] != "0000" }
}


def get_csv_4_args(channel) {
    return channel
        .map { taxid, csv, fastq1, fastq2 -> [csv, fastq1] }
        .transpose()
        .map { csv, fastq1 -> [fastq1.getSimpleName().replaceAll("_1", ""), csv] }
        .filter { it[0] != "0000" }
}


def get_busco_channel(sequences, busco_downloads, dataset_name, mode, priority) {
    return sequences.combine(busco_downloads).map {id, files, busco -> [id, files, dataset_name, mode, busco, priority]}
}


workflow DOWNLOAD_DATA {
    take:
        all
        all_without_models
        all_without_models_with_ncbi_key
        dataset_name

    main:
        // Download proteomes for outgroups
        DOWNLOAD_PROTEOMES(all)
        DOWNLOAD_REFSEQ_PROTEOMES(all)

        transpose_channel_3_args(DOWNLOAD_PROTEOMES.out)
            .concat(transpose_channel_3_args(DOWNLOAD_REFSEQ_PROTEOMES.out))
            .set { proteomes_unreformated } 
        
        REFORMAT_PROTEOMES(proteomes_unreformated)

        // Download genomes for the studies group 
        DOWNLOAD_GENOMES(all_without_models)

        // Download tsa for the studies group
        DOWNLOAD_TSA(all_without_models_with_ncbi_key)
        REFORMAT_TSA(transpose_channel_3_args(DOWNLOAD_TSA.out))
        CD_HIT_EST_TSA(REFORMAT_TSA.out)

        // Download SRA for the studies group
        DOWNLOAD_SRA(all_without_models_with_ncbi_key)
        // TRIMMOMATIC(transpose_channel_4_args(DOWNLOAD_SRA.out))
        TRINITY(transpose_channel_4_args(DOWNLOAD_SRA.out))
        REFORMAT_SRA(TRINITY.out)
        CD_HIT_EST_SRA(REFORMAT_SRA.out)

        // get the taxid for each set of sequences
        get_csv_3_args(DOWNLOAD_PROTEOMES.out)
            .concat(get_csv_3_args(DOWNLOAD_GENOMES.out), get_csv_3_args(DOWNLOAD_TSA.out), get_csv_4_args(DOWNLOAD_SRA.out))
            .set { csv }

        GET_TAXID(csv)

        GET_TAXID.out
            .map { id, taxid -> [id, file(taxid).readLines()[0]] }
            .set { taxid }

        // Download BUSCO datasets
        DOWNLOAD_BUSCO_DATASETS( Channel.from(dataset_name) )

        // concat all the single copy orthologs find for each specie
        GET_MARKER_NUMBER(DOWNLOAD_BUSCO_DATASETS.out.map { busco_dataset -> [dataset_name, busco_dataset]}) // get the number of marker in busco dataset

        // Combine proteomes, genomes, transcripts and tsa

        get_busco_channel(REFORMAT_PROTEOMES.out, DOWNLOAD_BUSCO_DATASETS.out, dataset_name, "proteins", 1)
            .concat(get_busco_channel(transpose_channel_3_args(DOWNLOAD_GENOMES.out), DOWNLOAD_BUSCO_DATASETS.out, dataset_name, "genome", 2))
            .concat(get_busco_channel(CD_HIT_EST_TSA.out, DOWNLOAD_BUSCO_DATASETS.out, dataset_name, "transcriptome", 3))
            .concat(get_busco_channel(CD_HIT_EST_SRA.out, DOWNLOAD_BUSCO_DATASETS.out, dataset_name, "transcriptome", 3))
            .set { busco_channel }

    emit:
        busco_channel   = busco_channel
        busco_datasets  = DOWNLOAD_BUSCO_DATASETS.out
        marker_number   = GET_MARKER_NUMBER.out.map { it[1].readLines()[0]}
        taxid           = taxid
}
