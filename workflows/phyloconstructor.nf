include { DOWNLOAD_PROTEOMES            } from "../modules/local/download_proteomes"
include { DOWNLOAD_GENOMES              } from "../modules/local/download_genomes"
include { DOWNLOAD_TSA                  } from "../modules/local/download_tsa"
include { DOWNLOAD_SRA                  } from "../modules/local/download_sra"
include { TRINITY                       } from "../modules/local/trinity"
include { BUSCO as BUSCO_PROTEINS       } from "../modules/local/busco"
include { BUSCO as BUSCO_GENOMES        } from "../modules/local/busco"
include { BUSCO as BUSCO_TRANSCRIPTOMES } from "../modules/local/busco"
include { DOWNLOAD_BUSCO_DATASETS       } from "../modules/local/download_busco_datasets"
include { REFORMAT_FASTA                } from "../modules/local/reformat_fasta"
include { GET_TAXID                     } from "../modules/local/get_taxid"


params.taxid = "61985"
params.outgroup_uniprot = "7227,6945,6669"
params.outgroup_all = "6850,438506,1519145"
params.threshold = 70
params.ncbi_api_key = "01679e044efe7ad60b87d93a1bb9085f7b09"
params.busco_dataset = "arthropoda_odb10"


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
//         .filter { it[0] != "0000" }
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


def get_busco_channel(sequences, busco_downloads, busco_dataset, mode) {
    return sequences.combine(busco_downloads).map {id, files, busco -> [id, files, busco_dataset, mode, busco]}
}

workflow PHYLOCONSTRUCTOR {
    model_outgroups = Channel.from(params.outgroup_uniprot.split(','))
    short_distance_outgroups = Channel.from(params.outgroup_all.split(','))
    taxid = Channel.from(params.taxid)

    all = taxid.concat(model_outgroups, short_distance_outgroups)
    all_without_models = taxid.concat(short_distance_outgroups)
    all_without_models_with_ncbi_key = all_without_models.map { taxid -> [taxid, params.ncbi_api_key] }

    // Download proteomes for outgroups
    DOWNLOAD_PROTEOMES(all)
    proteomes = transpose_channel_3_args(DOWNLOAD_PROTEOMES.out)

    // Download genomes for the studies group 
    DOWNLOAD_GENOMES(all_without_models)
    genomes = transpose_channel_3_args(DOWNLOAD_GENOMES.out)

    // Download tsa for the studies group
    DOWNLOAD_TSA(all_without_models_with_ncbi_key)
    tsa = transpose_channel_3_args(DOWNLOAD_TSA.out)
    REFORMAT_FASTA(tsa)

    // Download SRA for the studies group
    DOWNLOAD_SRA(all_without_models_with_ncbi_key)
    sra = transpose_channel_4_args(DOWNLOAD_SRA.out)
    // TRIMMOMATIC(sra)
    // TRINITY(TRIMMOMATIC.out)

    // // concatenate all transcripts files
    // REFORMAT_FASTA.out.concat(TRINITY.out).set { transcripts }

    // BUSCO
    DOWNLOAD_BUSCO_DATASETS( Channel.from(params.busco_dataset) )
    BUSCO_PROTEINS(get_busco_channel(proteomes, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "proteins"))
    BUSCO_GENOMES(get_busco_channel(genomes, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "genome"))
    // // BUSCO_TRANSCRIPTOMES(get_busco_channel(transcripts, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "transcriptome"))
    BUSCO_TRANSCRIPTOMES(get_busco_channel(REFORMAT_FASTA.out, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "transcriptome"))

    // // get the taxid for each file
    // get_csv_3_args(DOWNLOAD_PROTEOMES)
    //     .concat(get_csv_3_args(DOWNLOAD_GENOMES), get_csv_4_args(DOWNLOAD_TSA), get_csv_3_args(DOWNLOAD_SRA))
    //     .set { csv }
    // GET_TAXID(csv)

    // GET_TAXID.map { id, taxid -> id, file(taxid).readLines()[0] }.set { taxid }
    
    // // concat all the single copy orthologs find for each specie
    // GET_MARKER_NUMBER(DOWNLOAD_BUSCO_DATASETS.out) // get the number of marker in busco dataset

    // taxid.join(BUSCO_PROTEINS.out.busco_sequences)
    //     .concat(taxid.join(taxid.join(BUSCO_TRANSCRIPTOMES.out.busco_sequences, BUSCO_GENOMES.out.busco_sequences)))
    //     .map { id, taxid, busco -> [taxid, id, busco] }
    //     .groupTuple(0)
    //     .combine(GET_MARKER_NUMBER.out.readLines()[0])
    //     .set ( all_busco )

    // GET_SINGLE_COPY_ORTHOLOGS(all_busco)
    
    // // selection of species with complete percentage of single copy orthologs >= threshold
    // GET_SINGLE_COPY_ORTHOLOGS.out
    //     .map (taxid, sequences, csv -> [ taxid, sequences, csv.readLines()[0].split(',')[1].toInteger() ] )
    //     .filter ( it[2 >= 70 ])
    //     .map { taxid, sequences, csv -> [taxid, sequences] }
    //     .set { good_busco }
    
    // // concatenate all the single copy orthologs for the selected species
    // GROUP_BY_ORTHOLOGS(good_busco)
}

