include { FETCH_UNIPROT_PROTEOMES     } from '../modules/prepare_inputs/fetch_uniprot_proteomes'
include { DOWNLOAD_UNIPROT_PROTEOMES  } from '../modules/prepare_inputs/download_uniprot_proteomes'
include { FETCH_NCBI_ASSEMBLIES       } from '../modules/prepare_inputs/fetch_ncbi_assemblies'
include { DOWNLOAD_NCBI_ASSEMBLIES    } from '../modules/prepare_inputs/download_ncbi_assemblies'
include { GENERATE_NCBI_PROTEOMES     } from '../modules/prepare_inputs/generate_ncbi_proteomes'
include { FETCH_TSA_TRANSCRIPTOMES    } from '../modules/prepare_inputs/fetch_tsa_transcriptomes'
include { DOWNLOAD_TSA_TRANSCRIPTOMES } from '../modules/prepare_inputs/download_tsa_transcriptomes'
include { FETCH_SRA_READS             } from '../modules/prepare_inputs/fetch_sra_reads'
include { DOWNLOAD_SRA_READS          } from '../modules/prepare_inputs/download_sra_reads'
include { RNASPADES                   } from '../modules/prepare_inputs/rnaspades'
include { FILTER_TRANSCRIPTS          } from '../modules/prepare_inputs/filter_transcripts'
include { TRANSDECODER                } from '../modules/prepare_inputs/transdecoder'

workflow PREPARE_INPUTS {

    main:

    /*
     * Taxids from the focal group.
     * These taxids are used to query external databases.
     */
    Channel.from(params.group_taxid)
        .set { study_group_taxid }

    /*
     * User-provided outgroup proteomes.
     *
     * Expected CSV columns:
     *   taxid,name,fasta
     *
     * Species names are reformatted to lowercase identifiers
     * compatible with file names and downstream process tags.
     */
    Channel.fromPath(params.outgroups)
        .splitText()
        .filter { !it.startsWith('taxid') }
        .filter { it.trim() }
        .map { line ->
            def cols = line.split(',')

            def taxid = cols[0].trim()

            def specie = cols[1].trim()
                .replaceAll(/\s+/, '-')
                .replaceAll(/\./, '-')
                .toLowerCase()

            def file = cols[2].trim()

            [ taxid, specie, file ]
        }
        .set { model_outgroups }

    /*
     * User-provided proteomes from the focal group.
     *
     * Expected CSV columns:
     *   taxid,name,fasta
     */
    Channel.fromPath(params.group_species_csv)
        .splitText()
        .filter { !it.startsWith('taxid') }
        .filter { it.trim() }
        .map { line ->
            def cols = line.split(',')

            def taxid = cols[0].trim()

            def specie = cols[1]
                .trim()
                .replaceAll(/\s+/, '-')
                .replaceAll(/\./, '-')
                .toLowerCase()

            def file = cols[2].trim()

            [ taxid, specie, file ]
        }
        .set { model_group_species }

    /*
     * Taxids to query in external databases.
     *
     * This includes:
     *   - focal group taxids
     *   - outgroup taxids
     *
     * Duplicates are removed and values are sorted to make
     * downstream execution deterministic across runs.
     */
    study_group_taxid
        .concat(model_outgroups.map { it[0] })
        .distinct()
        .toSortedList()
        .flatMap { it }
        .set { taxids_to_fetch }

    /*
     * Combine all user-provided proteomes.
     *
     * Each tuple is formatted as:
     *   taxid, species, data_id, fasta
     *
     * The data_id is set to "user_proteomes" because these files
     * are directly provided by the user rather than downloaded.
     */
    model_outgroups
        .concat(model_group_species)
        .map { taxid, specie, fasta ->
            [ taxid, specie, 'user_proteomes', file(fasta) ]
        }
        .set { user_proteomes_raw }

    /*
     * Add BUSCO metadata to user-provided proteomes.
     *
     * Final tuple format:
     *   taxid, species, data_id, path, source, busco_mode, busco_dataset
     */
    user_proteomes = Utils.add_source_and_busco_dataset(
        user_proteomes_raw,
        'user_proteomes',
        'proteins',
        params.busco_lineage
    )

    /*
     * Optionally fetch and download UniProt proteomes.
     */
    if (params.download_uniprot == true) {

        /*
         * Fetch UniProt proteome IDs associated with each taxid.
         */
        FETCH_UNIPROT_PROTEOMES(taxids_to_fetch)

        /*
         * Parse the fetched ID files into tuples:
         *   taxid, species, proteome_id
         */
        proteome_ids = Utils.parse_ids_file(FETCH_UNIPROT_PROTEOMES.out.ids)

        /*
         * Download UniProt proteomes from their proteome IDs.
         */
        DOWNLOAD_UNIPROT_PROTEOMES(proteome_ids)

        /*
         * Add BUSCO metadata to downloaded UniProt proteomes.
         */
        uniprot_proteomes = Utils.add_source_and_busco_dataset(
            DOWNLOAD_UNIPROT_PROTEOMES.out,
            'uniprot_proteomes',
            'proteins',
            params.busco_lineage
        )

    } else {
        uniprot_proteomes = Channel.empty()
    }

    /*
     * Optionally fetch NCBI assemblies and/or NCBI proteomes.
     *
     * Assemblies are also fetched when SRA reads are requested, because
     * genome assemblies may be needed as references for downstream analyses.
     */
    if (params.download_ncbi_assemblies == true || params.download_ncbi_proteomes == true ) {

        /*
         * Fetch NCBI assembly and proteome IDs associated with each taxid.
         */
        FETCH_NCBI_ASSEMBLIES(taxids_to_fetch)

        /*
         * Parse NCBI assembly accessions.
         */
        ncbi_assembly_ids_raw = Utils.parse_ids_file(
            FETCH_NCBI_ASSEMBLIES.out.assembly_ids
        )

        /*
         * Prefer RefSeq assemblies when available.
         *
         * For each species:
         *   - keep GCF accessions if at least one exists
         *   - otherwise keep GCA accessions
         */
        ncbi_assembly_ids = Utils.keep_refseq_if_available(ncbi_assembly_ids_raw)

        /*
         * Parse NCBI proteome IDs.
         */
        ncbi_proteome_ids = Utils.parse_ids_file(
            FETCH_NCBI_ASSEMBLIES.out.proteomes_ids
        )

        /*
         * Download NCBI genome assemblies.
         */
        DOWNLOAD_NCBI_ASSEMBLIES(ncbi_assembly_ids)

        /*
        * Add BUSCO metadata to downloaded genome assemblies.
        *
        * DOWNLOAD_NCBI_ASSEMBLIES emits both genome FASTA and GFF.
        * BUSCO only needs the genome FASTA here, so the GFF is ignored.
        */
        ncbi_assemblies = Utils.add_source_and_busco_dataset(
            DOWNLOAD_NCBI_ASSEMBLIES.out.map { taxid, specie, accession, fna, gff -> [ taxid, specie, accession, fna ]},
            'ncbi_assemblies',
            'genome',
            params.busco_lineage
        )

        /*
         * Optionally download NCBI proteomes.
         */
        if (params.download_ncbi_proteomes == true) {

            DOWNLOAD_NCBI_ASSEMBLIES.out
                .filter { taxid, specie, accession, fna, gff -> gff.size() > 0 }
                .set { ncbi_assemblies_with_gff }

            GENERATE_NCBI_PROTEOMES(ncbi_assemblies_with_gff)    
            /*
             * Add BUSCO metadata to downloaded NCBI proteomes.
             */
            ncbi_proteomes = Utils.add_source_and_busco_dataset(
                GENERATE_NCBI_PROTEOMES.out.fasta,
                'ncbi_proteomes',
                'proteins',
                params.busco_lineage
            )

        } else {
            ncbi_proteomes = Channel.empty()
        }

    } else {
        ncbi_assemblies = Channel.empty()
        ncbi_proteomes = Channel.empty()
    }
    

        /*
     * Optionally fetch and download TSA transcriptomes.
     *
     * TSA transcriptomes are downloaded as transcript FASTA files,
     * then later filtered and translated into protein sequences.
     */
    if (params.download_tsa_transcriptomes == true) {

        FETCH_TSA_TRANSCRIPTOMES(taxids_to_fetch)

        tsa_ids = Utils.parse_ids_file(FETCH_TSA_TRANSCRIPTOMES.out.ids)

        DOWNLOAD_TSA_TRANSCRIPTOMES(tsa_ids)

        tsa_transcriptomes = Utils.add_source(
            DOWNLOAD_TSA_TRANSCRIPTOMES.out,
            'tsa_transcriptomes'
        )

    } else {
        tsa_transcriptomes = Channel.empty()
    }


    /*
     * Optionally fetch SRA reads and assemble transcriptomes with RNA-SPAdes.
     */
    if (params.download_sra_reads == true) {

        FETCH_SRA_READS(study_group_taxid)

        sra_ids = Utils.parse_ids_file(FETCH_SRA_READS.out.ids)

        DOWNLOAD_SRA_READS(sra_ids)

        DOWNLOAD_SRA_READS.out.sra_reads
            .map { taxid, specie, sra_id, reads1, reads2 ->
                [ taxid, specie, reads1, reads2 ]
            }
            .groupTuple(by: [0, 1])
            .set { reads_inputs }

        /*
         * TODO:
         * Add read subsampling if the number of reads is too high.
         */
        RNASPADES(reads_inputs)

        RNASPADES.out
            .map { taxid, specie, fasta ->
                [ taxid, specie, 'sra', fasta, 'sra_transcriptomes' ]
            }
            .set { rna_spades_transcriptomes }

    } else {
        rna_spades_transcriptomes = Channel.empty()
    }


    /*
     * Build transcript-derived proteomes.
     *
     * TSA transcriptomes and RNA-SPAdes transcriptomes are first filtered,
     * then translated into protein sequences with TransDecoder.
     */
    tsa_transcriptomes
        .concat(rna_spades_transcriptomes)
        .distinct { it[0..2].join('|') }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2]
        }
        .flatMap { it }
        .set { transcripts_inputs }

    FILTER_TRANSCRIPTS(transcripts_inputs)

    TRANSDECODER(FILTER_TRANSCRIPTS.out)

    TRANSDECODER.out
        .map { taxid, specie, data_id, path, source ->
            [ taxid, specie, data_id, path, source, 'proteins', params.busco_lineage ]
        }
        .distinct { it[0..2].join('|') }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2]
        }
        .flatMap { it }
        .set { all_transcripts }


    /*
     * Merge all available BUSCO inputs.
     *
     * Input sources currently included:
     *   - user-provided proteomes
     *   - UniProt proteomes
     *   - NCBI genome assemblies
     *   - NCBI proteomes
     *
     * Duplicate records are removed using:
     *   taxid + species + data_id
     *
     * Records are sorted to make BUSCO input order reproducible.
     */
    uniprot_proteomes
        .concat(user_proteomes, ncbi_assemblies, ncbi_proteomes, all_transcripts)
        .distinct { it[0..2].join('|') }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2] ?: a[4] <=> b[4]
        }
        .flatMap { it }
        .set { all_busco_inputs }

    emit:
    /*
     * Final BUSCO input channel.
     *
     * Tuple format:
     *   taxid, species, data_id, path, source, busco_mode, busco_dataset
     */
    all_busco_inputs = all_busco_inputs
}