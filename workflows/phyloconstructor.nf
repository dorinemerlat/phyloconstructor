include { DOWNLOAD_PROTEOMES            } from "../modules/local/download_proteomes"
include { DOWNLOAD_GENOMES              } from "../modules/local/download_genomes"
include { DOWNLOAD_TSA                  } from "../modules/local/download_tsa"
include { DOWNLOAD_SRA                  } from "../modules/local/download_sra"
include { ASSEMBLE_SRA                  } from "../modules/local/assemble_sra"
include { BUSCO as BUSCO_PROTEINS       } from "../modules/local/busco"
include { BUSCO as BUSCO_GENOMES        } from "../modules/local/busco"
include { BUSCO as BUSCO_TRANSCRIPTOMES } from "../modules/local/busco"
include { DOWNLOAD_BUSCO_DATASETS       } from "../modules/local/download_busco_datasets"
include { REFORMAT_FASTA                } from "../modules/local/reformat_fasta"

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
        .map { taxid, fasta -> [fasta.getBaseName(), fasta] }
}

def transpose_channel_4_args(channel) {
    return channel
        .map { taxid, csv, fastq1, fastq2 -> [taxid, fastq1, fastq2] }
        .transpose()
        .map { taxid, fastq1, fastq2 -> [fasta.getBaseName(), fastq1, fastq2] }
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
    TRIMMOMATIC(sra)
    TRINITY(TRIMMOMATIC.out)

    // concatenate all transcripts files
    REFORMAT_FASTA.out.concat(TRINITY.out).sey { transcripts }

    // BUSCO
    DOWNLOAD_BUSCO_DATASETS( Channel.from(params.busco_dataset) )

    BUSCO_PROTEINS(get_busco_channel(proteomes, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "proteins"))
    BUSCO_GENOMES(get_busco_channel(genomes, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "genome"))
    BUSCO_TRANSCRIPTOMES(get_busco_channel(transcripts, DOWNLOAD_BUSCO_DATASETS.out, params.busco_dataset, "transcriptome"))

    // Add more processes as needed for alignment, concatenation, etc.

}

