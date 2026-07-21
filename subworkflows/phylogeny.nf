include { EXTRACT_BUSCO_FASTA         } from '../modules/phylogeny/extract_busco_fasta'
include { MAFFT                       } from '../modules/phylogeny/mafft'
include { TRIMAL                      } from '../modules/phylogeny/trimal'
include { PHYKIT_CONCAT               } from '../modules/phylogeny/phykit_concat'
include { IQTREE_SUPERMATRIX          } from '../modules/phylogeny/iqtree_supermatrix'
include { IQTREE_INDIVIDUAL_GENE_TREE } from '../modules/phylogeny/iqtree_individual_gene_tree'
include { ASTRAL                      } from '../modules/phylogeny/astral'


workflow PHYLOGENY {

    take:
    selected_orthogroups
    single_copy_busco_sequences
    multi_copy_busco_sequences

    main:

    /*
     * Convert the user-defined threshold combinations to integers
     * and store them as a Set for efficient membership testing.
     *
     * Parameter validation is performed once in main.nf before
     * the pipeline starts.
     */
    def thresholds_to_test = params.phylogeny_thresholds
        .collect { threshold ->
            [
                threshold[0] as Integer,
                threshold[1] as Integer
            ]
        }
        .toSet()


    /*
     * Associate each selected orthogroup table with:
     *
     *   - its dataset label;
     *   - the BUSCO completeness threshold;
     *   - the gene occupancy threshold;
     *   - the single-copy BUSCO sequences;
     *   - the multi-copy BUSCO sequences.
     */
    selected_orthogroups
        .flatMap { label, tables ->
            tables.collect { table ->
                [label, table]
            }
        }
        .map { label, table ->

            def match = (
                table.name =~
                /orthogroups_busco-c_(\d+)_gene-occupancy_(\d+)\.tsv/
            )

            if (!match.matches()) {
                error "Cannot parse thresholds from file name: ${table}"
            }

            def busco_c_threshold =
                match[0][1] as Integer

            def gene_occupancy_threshold =
                match[0][2] as Integer

            [
                label,
                busco_c_threshold,
                gene_occupancy_threshold,
                table
            ]
        }
        .filter {
            label,
            busco_c_threshold,
            gene_occupancy_threshold,
            table ->

            [
                busco_c_threshold,
                gene_occupancy_threshold
            ] in thresholds_to_test
        }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?:
            a[1] <=> b[1] ?:
            a[2] <=> b[2]
        }
        .flatMap { it }
        .combine(
            single_copy_busco_sequences.map {
                [it]
            }
        )
        .combine(
            multi_copy_busco_sequences.map {
                [it]
            }
        )
        .map {
            label,
            busco_c_threshold,
            gene_occupancy_threshold,
            table,
            single_copy_sequences,
            multi_copy_sequences ->

            def job_name = [
                "busco-c_${busco_c_threshold}",
                "gene-occupancy_${gene_occupancy_threshold}"
            ].join('_')

            [
                label,
                job_name,
                table,
                single_copy_sequences,
                multi_copy_sequences
            ]
        }
        .set { orthogroups_with_busco_sequences }


    /*
     * Extract one FASTA file per retained BUSCO orthogroup.
     */
    EXTRACT_BUSCO_FASTA(
        orthogroups_with_busco_sequences
    )


    /*
     * Convert the FASTA lists into one channel item per orthogroup.
     *
     * Expected tuple:
     *
     *   [label, job_name, orthogroup, fasta]
     */
    EXTRACT_BUSCO_FASTA.out
        .flatMap { label, job_name, fasta_list ->

            fasta_list.collect { fasta ->

                def orthogroup =
                    fasta.baseName.tokenize('_')[-1]

                [
                    label,
                    job_name,
                    orthogroup,
                    fasta
                ]
            }
        }
        .distinct {
            it[0..2].join('|')
        }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?:
            a[1] <=> b[1] ?:
            a[2] <=> b[2]
        }
        .flatMap { it }
        .set { all_orthogroup_fastas }


    /*
     * Align each orthogroup independently with MAFFT.
     */
    MAFFT(
        all_orthogroup_fastas
    )


    /*
     * Trim each multiple sequence alignment with trimAl.
     */
    TRIMAL(
        MAFFT.out
    )


    /*
     * Group trimmed alignments by dataset and threshold combination,
     * then sort them by orthogroup before concatenation.
     *
     * Output tuple:
     *
     *   [label, job_name, list_of_alignments]
     */
    TRIMAL.out
        .groupTuple(by: [0, 1])
        .map {
            label,
            job_name,
            orthogroups,
            alignments ->

            def sorted_alignments =
                [orthogroups, alignments]
                    .transpose()
                    .sort { it[0] }

            [
                "trimal_${label}",
                job_name,
                sorted_alignments.collect {
                    it[1]
                }
            ]
        }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?:
            a[1] <=> b[1]
        }
        .flatMap { it }
        .set { all_trimal_alignments }


    /*
     * Concatenate the trimmed alignments into a supermatrix and generate
     * the corresponding partition file.
     */
    PHYKIT_CONCAT(
        all_trimal_alignments
    )


    /*
     * Infer the concatenation-based maximum-likelihood tree.
     */
    IQTREE_SUPERMATRIX(
        PHYKIT_CONCAT.out.main_output
    )


    /*
     * Prepare trimmed alignments for individual gene-tree inference.
     *
     * Input tuple:
     *
     *   [label, job_name, orthogroup, alignment]
     *
     * The label is prefixed with "trimal_" to identify the trimming
     * method in downstream output names.
     */
    TRIMAL.out
        .map {
            label,
            job_name,
            orthogroup,
            alignment ->

            [
                "trimal_${label}",
                job_name,
                orthogroup,
                alignment
            ]
        }
        .set { trimmed_alignments_for_gene_trees }


    /*
     * Infer one maximum-likelihood tree per orthogroup.
     */
    IQTREE_INDIVIDUAL_GENE_TREE(
        trimmed_alignments_for_gene_trees
    )


    /*
     * Group individual gene trees by dataset and threshold combination
     * before coalescent-based species-tree inference.
     */
    IQTREE_INDIVIDUAL_GENE_TREE.out.treefile
        .groupTuple(by: [0, 1])
        .set { gene_trees_by_label }


    /*
     * Infer the coalescent-based species tree with ASTRAL.
     */
    ASTRAL(
        gene_trees_by_label
    )
}