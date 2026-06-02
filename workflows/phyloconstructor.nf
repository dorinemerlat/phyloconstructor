include { FETCH_UNIPROT_PROTEOMES       }  from '../modules/fetch_uniprot_proteomes'
include { DOWNLOAD_UNIPROT_PROTEOMES    }  from '../modules/download_uniprot_proteomes'
include { FETCH_NCBI_ASSEMBLIES         }  from '../modules/fetch_ncbi_assemblies'
include { DOWNLOAD_NCBI_ASSEMBLIES      }  from '../modules/download_ncbi_assemblies'
include { DOWNLOAD_NCBI_PROTEOMES       }  from '../modules/download_ncbi_proteomes'
include { FETCH_TSA_TRANSCRIPTOMES      }  from '../modules/fetch_tsa_transcriptomes'
include { DOWNLOAD_TSA_TRANSCRIPTOMES   }  from '../modules/download_tsa_transcriptomes'
include { FETCH_SRA_READS               }  from '../modules/fetch_sra_reads'
include { DOWNLOAD_SRA_READS            }  from '../modules/download_sra_reads'
include { STAR                          }  from '../modules/star'
include { STRINGTIE                     }  from '../modules/stringtie'
include { RNASPADES                     }  from '../modules/rnaspades'
include { FILTER_TRANSCRIPTS            }  from '../modules/filter_transcripts'
include { TRANSDECODER                  }  from '../modules/transdecoder'
include { BUSCO                         }  from '../modules/busco'
include { FILTER_COMPLETE_BUSCO_HITS    }  from '../modules/filter_complete_busco_hits'
include { SELECT_BEST_BUSCO_HITS        }  from '../modules/select_best_busco_hits'
include { SELECT_ORTHOGROUPS            }  from '../modules/select_orthogroups'
include { EXTRACT_BUSCO_FASTA           }  from '../modules/extract_busco_fasta'
include { MAFFT                         }  from '../modules/mafft'
include { TRIMAL                        }  from '../modules/trimal'
include { AMAS                          }  from '../modules/amas'
include { IQTREE_SUPERMATRIX            }  from '../modules/iqtree_supermatrix'
include { IQTREE_INDIVIDUAL_GENE_TREE   }  from '../modules/iqtree_individual_gene_tree'
include { ASTRAL                        }  from '../modules/astral'

def parse_ids_file(ch) {
    ch.flatMap { taxid, ids_file ->

        ids_file.text
            .readLines()
            .drop(1)
            .findAll { it.trim() }
            .collect { line ->

                def cols = line.split('\t')

                if (cols.size() < 3) {
                    log.warn "Skipping malformed line in ${ids_file}: ${line}"
                    return null
                }

                def organism_taxid = cols[0].trim()
                def organism_name = cols[1].trim().replaceAll(/\./, '-')
                def data_id = cols[2].trim()

                tuple(organism_taxid, organism_name, data_id)
            }
        }
        .distinct()
}

def keep_refseq_if_available(ch) {
    ch
        .groupTuple(by: [0, 1])
        .map { taxid, organism_name, data_ids ->

            def ids = data_ids.collect { it.toString().trim() }.unique()

            def gcf_ids = ids.findAll { it.startsWith("GCF_") }
            def gca_ids = ids.findAll { it.startsWith("GCA_") }

            def selected_ids

            if (gcf_ids) {
                selected_ids = gcf_ids
            } else {
                selected_ids = gca_ids
            }

            selected_ids.collect { data_id ->
                tuple(taxid, organism_name, data_id)
            }
        }
        .flatMap { it }
}

def add_source_and_busco_dataset(ch, source, mode, dataset) {
    ch.map { taxid, specie, id, path -> [taxid, specie, id, path, source, mode, dataset] }
}

def add_source(ch, source) {
    ch.map { taxid, specie, id, path -> [taxid, specie, id, path, source] }
}
workflow PHYLOCONSTRUCTOR {

    main:
    // get all taxids
    // Studied taxid group
    Channel.from(params.taxid)
        .set { taxid }

    // Outgroups 
    Channel
        .fromPath(params.outgroups)
        .splitText()
        .filter { !it.startsWith('taxid') }
        .map { line -> 
            cols = line.split(',')
            taxid = cols[0].trim()
            def specie = cols[1].trim()
                .replaceAll(/\s+/, '-')
                .toLowerCase()
            def file = cols[2].trim()
            [taxid, specie, file]}
        .set { model_outgroups }

    // Interest group 
    Channel
        .fromPath(params.group_species_csv)
        .splitText()
        .filter { !it.startsWith('taxid') }
        .map { line -> 
            cols = line.split(',')
            taxid = cols[0].trim()
            def specie = cols[1].trim()
            specie = specie.replaceAll(/\s+/, '-')
            specie = specie.toLowerCase()
            def file = cols[2].trim()
            [taxid, specie, file]}
        .set { model_group_species }

    // Concat studied taxid group and outgroups (remove duplicates)
    taxid
        .concat(model_outgroups.map { it[0] }, model_group_species.map { it[0] })
        .distinct()
        .set { all_taxids }

    // Concat studied taxid group and outgroups with their corresponding files
    model_outgroups.concat(model_group_species)
        .filter { it[2] } // filter for lines where file is not empty
        .map { taxid, specie, fasta -> [taxid, specie, 'user_proteomes', file(fasta)] }
        .set { user_proteomes }

    user_proteomes = add_source_and_busco_dataset(user_proteomes, 'user_proteomes', 'proteins', params.busco_lineage)

    if (params.download_uniprot == true) {
        // fetch for proteomes in uniprot for the taxid of interest and the outgroups and download the corresponding proteomes
        FETCH_UNIPROT_PROTEOMES(all_taxids)
        proteome_ids = parse_ids_file(FETCH_UNIPROT_PROTEOMES.out.ids)
        DOWNLOAD_UNIPROT_PROTEOMES(proteome_ids)
        uniprot_proteomes = add_source_and_busco_dataset(DOWNLOAD_UNIPROT_PROTEOMES.out, 'uniprot_proteomes', 'proteins', params.busco_lineage)
    } else {
        uniprot_proteomes = Channel.empty()
    }

    if (params.download_ncbi_assemblies == true || params.download_ncbi_proteomes == true || params.download_sra_reads == true) {
        FETCH_NCBI_ASSEMBLIES(all_taxids)

        ncbi_assembly_ids_raw = parse_ids_file(FETCH_NCBI_ASSEMBLIES.out.assembly_ids)
        ncbi_assembly_ids = keep_refseq_if_available(ncbi_assembly_ids_raw)

        ncbi_proteome_ids = parse_ids_file(FETCH_NCBI_ASSEMBLIES.out.proteomes_ids)

        if (params.download_ncbi_assemblies == true || params.download_sra_reads == true) {
            // fetch for genomes in NCBI for the taxid of interest and the outgroups and download the corresponding genomes
            DOWNLOAD_NCBI_ASSEMBLIES(ncbi_assembly_ids)
            ncbi_assemblies = add_source_and_busco_dataset(DOWNLOAD_NCBI_ASSEMBLIES.out, 'ncbi_assemblies', 'genome', params.busco_lineage)
        } else {
            ncbi_assemblies = Channel.empty()
        }
        
        if (params.download_ncbi_proteomes == true ) {
            // fetch for proteomes in NCBI for the taxid of interest and the outgroups and download the corresponding proteomes
            DOWNLOAD_NCBI_PROTEOMES(ncbi_proteome_ids)
            ncbi_proteomes = add_source_and_busco_dataset(DOWNLOAD_NCBI_PROTEOMES.out, 'ncbi_proteomes', 'proteins', params.busco_lineage)

        } else {
            ncbi_proteomes = Channel.empty()
        }
    }
 
    if (params.download_tsa_transcriptomes == true) {
        // fetch for whole asssembly assembly transcriptomes for the taxid of interest and the outgroups and download the corresponding transcriptomes
        FETCH_TSA_TRANSCRIPTOMES(all_taxids)
        tsa_ids = parse_ids_file(FETCH_TSA_TRANSCRIPTOMES.out.ids)
        DOWNLOAD_TSA_TRANSCRIPTOMES(tsa_ids)
        tsa_transcriptomes = add_source(DOWNLOAD_TSA_TRANSCRIPTOMES.out, 'tsa_transcriptomes')
    } else {
        tsa_transcriptomes = Channel.empty()
    }

    if (params.download_sra_reads == true) {
        // fetch for SRA reads for the taxid of interest and the outgroups and download the corresponding reads
        FETCH_SRA_READS(taxid)
        sra_ids = parse_ids_file(FETCH_SRA_READS.out.ids)
        DOWNLOAD_SRA_READS(sra_ids)
        DOWNLOAD_SRA_READS.out.sra_reads
            .map { taxid, specie, sra_id, reads1, reads2 -> [taxid, specie, reads1, reads2] }
            .groupTuple(by: [0, 1])
            .set { reads_inputs }
        
        // TODO: Ajouter un sous-echantillonage a si threshold de reads > 50M de reads pour eviter plantage et target a 50M de reads

        RNASPADES(reads_inputs)
        RNASPADES.out
            .map{taxid, specie, fasta -> [taxid, specie, 'sra', fasta, 'sra_transcriptomes']}
            .set { rna_spades_transcriptomes }
        
    } else {
        rna_spades_transcriptomes = Channel.empty()
    }

    // process transcripts 
    tsa_transcriptomes
        .concat(rna_spades_transcriptomes)
        .set { transcripts_inputs }
    
    FILTER_TRANSCRIPTS(transcripts_inputs)
    TRANSDECODER(FILTER_TRANSCRIPTS.out)
    TRANSDECODER.out
        .map { taxid, specie, data_id, path, source -> [taxid, specie, data_id, path, source, 'proteins', params.busco_lineage] }
        .set { all_transcripts }

    // merge all inputs for busco
    uniprot_proteomes
        .concat(user_proteomes, ncbi_assemblies, ncbi_proteomes, all_transcripts)
        .set { all_busco_inputs }

    BUSCO(all_busco_inputs)
    FILTER_COMPLETE_BUSCO_HITS(BUSCO.out.full_table) 
    FILTER_COMPLETE_BUSCO_HITS.out
        .map { taxid, specie, data_id, table, source -> [taxid, specie, table] }
        .groupTuple(by: 0)
        .map { taxid, species, tables -> [taxid, species[0], tables.unique()] }
        .set { all_complete_tables }
    SELECT_BEST_BUSCO_HITS(all_complete_tables)

    SELECT_BEST_BUSCO_HITS.out
        .map { taxid, specie, table -> ['all', table] }
        .groupTuple(by: 0)
        .map { id, tables -> [tables, params.busco_dataset_size] }
        .distinct()
        .set { all_busco_results }

    SELECT_ORTHOGROUPS(all_busco_results)

    def thresholds_to_test = [
        [70,70],
        [70,80],
        [80,80],
        [80,90],
        [70,90]
    ] as Set

    BUSCO.out.single_copy_busco_sequences
        .concat(BUSCO.out.multi_copy_busco_sequences)
        .map { taxid, specie, data_id, fasta, source -> fasta }
        .flatten()
        .unique()
        .collect()
        .map { [it] }
        .set { all_busco_sequences }

    SELECT_ORTHOGROUPS.out.tables
        .flatten()
        .map { file ->
            def m = (file.name =~ /orthogroups_busco-c_(\d+)_gene-occupancy_(\d+)\.tsv/)

            if (!m.matches()) {
                error "Cannot parse thresholds from file name: ${file}"
            }

            def busco_c_threshold = m[0][1] as Integer
            def gene_occupancy_threshold = m[0][2] as Integer

            tuple(busco_c_threshold, gene_occupancy_threshold, file)
        }
        .filter { busco_c_threshold, gene_occupancy_threshold, file ->
            [busco_c_threshold, gene_occupancy_threshold] in thresholds_to_test
        }
        // .combine(all_busco_sequences)
        // .map { busco_threshold, gene_occupancy_threshold, table, busco_sequences ->
        //     tuple(
        //         busco_threshold,
        //         gene_occupancy_threshold,
        //         table,
        //         busco_sequences.unique()
        //     )
        // }
        .map {busco_c_threshold, gene_occupancy_threshold, file -> [ "busco-c_${busco_c_threshold}_gene-occupancy_${gene_occupancy_threshold}", file] }
        .set { orthogroups_with_busco_sequences }

    EXTRACT_BUSCO_FASTA(orthogroups_with_busco_sequences)

    EXTRACT_BUSCO_FASTA.out
        .flatMap { job_name, fasta_list ->
            fasta_list.collect { fasta ->
                def orthogroup = fasta.baseName.split('_')[-1]
                [ job_name, orthogroup, fasta ]
            }
        }
        .set { all_orthogroup_fastas }

    MAFFT(all_orthogroup_fastas)
    TRIMAL(MAFFT.out)

    TRIMAL.out
        .map { job_name, orthogroup, aln -> [job_name, aln] }
        .groupTuple(by: [0])
        .set { all_alns }

    AMAS(all_alns)

    AMAS.out
        .combine(model_outgroups.map { it[1] }.collect().map { [it] })
        .map { job_name, aln, partition, model_outgroups -> [job_name, aln, partition, model_outgroups.join(',')] }
        .set { supermatrix_inputs }


    TRIMAL.out.set { individual_gene_tree_inputs }
    

    IQTREE_SUPERMATRIX(supermatrix_inputs)
    IQTREE_INDIVIDUAL_GENE_TREE(individual_gene_tree_inputs)

    ASTRAL(IQTREE_INDIVIDUAL_GENE_TREE.out.treefile.groupTuple(by: 0))
    // emit:
    // busco_results = BUSCO.out.results
}