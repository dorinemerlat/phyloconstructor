# Phyloconstructor

Phyloconstructor is a Nextflow pipeline for assembling phylogenomic datasets from heterogeneous sequence sources (NCBI, UniProt, SRA, assemblies and transcriptomes), extracting BUSCO markers, selecting orthogroups, aligning and trimming gene sets, and producing gene trees and supermatrix/summary species trees.

Key features
- Data acquisition: fetch proteomes, assemblies, transcriptomes and SRA reads.
- BUSCO-based marker extraction and reformatting.
- Orthogroup selection, clustering and filtering.
- Alignment (MAFFT), trimming (trimAl), gene-tree inference (IQ-TREE) and supermatrix building.
- Support for downstream species-tree inference (ASTRAL) and summary outputs.

Requirements
- Nextflow (DSL2) and Java
- Python 3
- Common bioinformatics tools: BUSCO, MAFFT, IQ-TREE, trimAl, CD-HIT, TransDecoder, rnaspades, STAR, StringTie

Quickstart
- Configure sample lists in `group_species.csv` and `outgroups.csv`, or set `params.*` in `nextflow.config`.
- Run the pipeline from the repository root:

```
module load nextflow singularity
nextflow run main.nf
```

On a SLURM cluster, you can also submit [nextflow.sbatch](nextflow.sbatch) instead of running the command manually.

Configuration & layout
- Main config: `nextflow.config` (project params and included module config files).
- Workflows: [workflows/phyloconstructor.nf](workflows/phyloconstructor.nf)
- Modules: see the `modules/` directory for individual processes (MAFFT, IQ-TREE, BUSCO, etc.).
- Input data: `data/group_species/` and `data/outgroups/` as referenced in the CSV files.

Outputs & logs
- Nextflow working dirs and process outputs are created under the workspace (`work/`, `cache/`, and `log_phyloconstructor/`).

Contributing & citation
- For issues or contributions, open a GitHub issue or pull request.
- If you use Phyloconstructor in a publication, please cite the repository and include relevant tool citations for core components (BUSCO, MAFFT, IQ-TREE).

License
- See repository for license information.

