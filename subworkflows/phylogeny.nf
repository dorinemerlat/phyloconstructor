include { EXTRACT_BUSCO_FASTA         } from '../modules/phylogeny/extract_busco_fasta'
include { MAFFT                       } from '../modules/phylogeny/mafft'
include { TRIMAL                      } from '../modules/phylogeny/trimal'
include { CLIPKIT                     } from '../modules/phylogeny/clipkit'
include { AMAS                        } from '../modules/phylogeny/amas'
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

    def thresholds_to_test = [
        [60, 80],
        [60, 90],
        [70, 80],
        [70, 90],
        [80, 80],
        [80, 90],
    ] as Set

    selected_orthogroups
        .flatMap { label, tables ->
            tables.collect { table ->
                [ label, table ]
            }
        }
        .map { label, file ->
            def m = (file.name =~ /orthogroups_busco-c_(\d+)_gene-occupancy_(\d+)\.tsv/)

            if (!m.matches()) {
                error "Cannot parse thresholds from file name: ${file}"
            }

            def busco_c_threshold = m[0][1] as Integer
            def gene_occupancy_threshold = m[0][2] as Integer

            [ label, busco_c_threshold, gene_occupancy_threshold, file ]
        }
        .filter { label, busco_c_threshold, gene_occupancy_threshold, file ->
            [ busco_c_threshold, gene_occupancy_threshold ] in thresholds_to_test
        }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2]
        }
        .flatMap { it }
        .combine(single_copy_busco_sequences.map { [it] })
        .combine(multi_copy_busco_sequences.map { [it] })
        .map { label, busco_c_threshold, gene_occupancy_threshold, table,
              single_copy_sequences, multi_copy_sequences ->

            def job_name = "busco-c_${busco_c_threshold}_gene-occupancy_${gene_occupancy_threshold}"

            [ label, job_name, table, single_copy_sequences, multi_copy_sequences ]
        }
        .set { orthogroups_with_busco_sequences }

    EXTRACT_BUSCO_FASTA(orthogroups_with_busco_sequences)

    EXTRACT_BUSCO_FASTA.out
        .flatMap { label, job_name, fasta_list ->
            fasta_list.collect { fasta ->
                def orthogroup = fasta.baseName.tokenize('_')[-1]
                [ label, job_name, orthogroup, fasta ]
            }
        }
        .distinct { it[0..2].join('|') }
        .toSortedList { a, b ->
            a[0] <=> b[0] ?: a[1] <=> b[1] ?: a[2] <=> b[2]
        }
        .flatMap { it }
        .set { all_orthogroup_fastas }

    MAFFT(all_orthogroup_fastas)

    TRIMAL(MAFFT.out)
    CLIPKIT(MAFFT.out)

    TRIMAL.out
        .groupTuple(by: [0, 1])
        .map { label, job_name, orthogroups, alns ->
            def sorted = [orthogroups, alns].transpose().sort { it[0] }
            [ 'trimal_' + label, job_name, sorted.collect { it[1] } ]
        }
        .toSortedList { a, b ->
            a[0] <=> b[0]
        }
        .flatMap { it }
        .set { all_alns_trimal }

    CLIPKIT.out
        .groupTuple(by: [0, 1])
        .map { label, job_name, orthogroups, alns ->
            def sorted = [orthogroups, alns].transpose().sort { it[0] }
            [ 'clipseq_' + label, job_name, sorted.collect { it[1] } ]
        }
        .toSortedList { a, b ->
            a[0] <=> b[0]
        }
        .flatMap { it }
        .set { all_alns_clipkit }

    all_alns_trimal.concat(all_alns_clipkit).set{all_alns} 

    PHYKIT_CONCAT(all_alns)

    IQTREE_SUPERMATRIX(PHYKIT_CONCAT.out.main_output)

    TRIMAL.out
        .map { label, job_name, orthogroups, alns -> ['trimal_' + label, job_name, orthogroups, alns ] }
        .set {trimal_out}
    CLIPKIT.out
        .map { label, job_name, orthogroups, alns -> ['clipseq_' + label, job_name, orthogroups, alns ] }
        .set {clipseq_out}
    trimal_out.concat(clipseq_out).set{clean_alns}
    IQTREE_INDIVIDUAL_GENE_TREE(clean_alns)

    IQTREE_INDIVIDUAL_GENE_TREE.out.treefile
        .groupTuple(by: [0, 1])
        .set { gene_trees_by_label }
    
    gene_trees_by_label.view()

    ASTRAL(gene_trees_by_label)
}
