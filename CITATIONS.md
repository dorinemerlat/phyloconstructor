# PhyloConstructor: Citations

When using PhyloConstructor, please cite the PhyloConstructor GitHub repository together with the publications associated with the software tools used in your analysis.

The exact tools executed depend on the selected input sources and workflow parameters.

## PhyloConstructor

GitHub repository:
https://github.com/dorinemerlat/phyloconstructor

## Workflow engine

### [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

> Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017;35(4):316–319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

### nf-metro

The workflow overview figure included in this repository was initially generated using **nf-metro** and subsequently refined manually.

Project:
https://github.com/nextflow-io/nf-metro

## Public sequence resources

PhyloConstructor can retrieve genomic, proteomic and transcriptomic data from UniProt and NCBI resources. Please cite the databases from which data were obtained.

### [UniProt](https://pubmed.ncbi.nlm.nih.gov/36408920/)

> The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res. 2023;51(D1):D523–D531. doi: 10.1093/nar/gkac1052. PubMed PMID: 36408920.

### [NCBI Datasets](https://doi.org/10.1038/s41597-024-03571-y)

> O'Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data. 2024;11:732. doi:10.1038/s41597-024-03571-y.

### [NCBI Sequence Read Archive](https://pubmed.ncbi.nlm.nih.gov/21062823/)

> Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. Nucleic Acids Res. 2011;39(Database issue):D19–D21. doi: 10.1093/nar/gkq1019. PubMed PMID: 21062823; PubMed Central PMCID: PMC3013647.

### [NCBI Transcriptome Shotgun Assembly database](https://www.ncbi.nlm.nih.gov/genbank/tsa/)

> National Center for Biotechnology Information. Transcriptome Shotgun Assembly Sequence Database. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information. Available from: https://www.ncbi.nlm.nih.gov/genbank/tsa/.

## Dataset completeness and orthologue selection

### [BUSCO](https://pubmed.ncbi.nlm.nih.gov/34320186/)

> Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Mol Biol Evol. 2021;38(10):4647–4654. doi: 10.1093/molbev/msab199. PubMed PMID: 34320186; PubMed Central PMCID: PMC8476166.

When reporting BUSCO results, the lineage dataset and its version should also be stated.

Example:

> BUSCO v6 was run using the `arthropoda_odb10` lineage dataset.

## Transcriptome reconstruction and protein prediction

These tools are used when PhyloConstructor constructs protein datasets from TSA transcriptomes or SRA RNA-seq reads.

### [SRA Toolkit](https://github.com/ncbi/sra-tools)

> National Center for Biotechnology Information. SRA Toolkit. Available from: https://github.com/ncbi/sra-tools.

### [BBMap](https://sourceforge.net/projects/bbmap/)

> Bushnell B. BBMap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Laboratory. Available from: https://sourceforge.net/projects/bbmap/.

### [SPAdes and RNA-SPAdes](https://pubmed.ncbi.nlm.nih.gov/31494669/)

> Bushmanova E, Antipov D, Lapidus A, Prjibelski AD. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience. 2019;8(9):giz100. doi: 10.1093/gigascience/giz100. PubMed PMID: 31510679; PubMed Central PMCID: PMC6751136.

The original SPAdes publication may also be cited:

> Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD, Pyshkin AV, Sirotkin AV, Vyahhi N, Tesler G, Alekseyev MA, Pevzner PA. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol. 2012;19(5):455–477. doi: 10.1089/cmb.2012.0021. PubMed PMID: 22506599; PubMed Central PMCID: PMC3342519.

### [TransDecoder](https://github.com/TransDecoder/TransDecoder)

> Haas BJ, Papanicolaou A. TransDecoder: find coding regions within transcripts. Available from: https://github.com/TransDecoder/TransDecoder.

## Multiple-sequence alignment and trimming

### [MAFFT](https://academic.oup.com/bioinformatics/article/32/13/1933/1743504)

> Kuraku S, Zmasek CM, Nishimura O, Katoh K. aLeaves facilitates on-demand exploration of metazoan gene family trees on MAFFT sequence alignment server with enhanced interactivity. Nucleic Acids Research. 2013;41(W1):W22-W28. doi:10.1093/nar/gkt389.

### [trimAl](https://pubmed.ncbi.nlm.nih.gov/19505945/)

> Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Bioinformatics. 2009;25(15):1972–1973. doi: 10.1093/bioinformatics/btp348. PubMed PMID: 19505945; PubMed Central PMCID: PMC2712344.

## Alignment concatenation and sequence utilities

### [PhyKIT](https://pubmed.ncbi.nlm.nih.gov/33560364/)

> Steenwyk JL, Buida TJ III, Li Y, Shen XX, Rokas A. PhyKIT: a broadly applicable UNIX shell toolkit for processing and analyzing phylogenomic data. Bioinformatics. 2021;37(16):2325–2331. doi: 10.1093/bioinformatics/btab096. PubMed PMID: 33892497; PubMed Central PMCID: PMC8370967.

### [SeqKit](https://pubmed.ncbi.nlm.nih.gov/27706213/)

> Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS One. 2016;11(10):e0163962. doi: 10.1371/journal.pone.0163962. PubMed PMID: 27706213; PubMed Central PMCID: PMC5051824.

### [CD-HIT](https://pubmed.ncbi.nlm.nih.gov/16731699/)

> Li W, Godzik A. Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics. 2006;22(13):1658–1659. doi: 10.1093/bioinformatics/btl158. PubMed PMID: 16731699.

### [AGAT](https://doi.org/10.5281/zenodo.3552717)

> Dainat J, Hereñú D, LucileSol, Pascal-Git. NBISweden/AGAT: AGAT. Zenodo. doi: 10.5281/zenodo.3552717.

## Phylogenetic inference

### [IQ-TREE](https://academic.oup.com/mbe/article/43/5/msag117/8669857)

> Minh BQ, Schmidt HA, Schrempf D, et al. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution. 2025.

### [ModelFinder](https://pubmed.ncbi.nlm.nih.gov/28481363/)

> Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. ModelFinder: fast model selection for accurate phylogenetic estimates. Nat Methods. 2017;14(6):587–589. doi: 10.1038/nmeth.4285. PubMed PMID: 28481363; PubMed Central PMCID: PMC5453245.

### [Ultrafast bootstrap](https://pubmed.ncbi.nlm.nih.gov/29077904/)

> Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. UFBoot2: improving the ultrafast bootstrap approximation. Mol Biol Evol. 2018;35(2):518–522. doi: 10.1093/molbev/msx281. PubMed PMID: 29077904; PubMed Central PMCID: PMC5850222.

## Coalescent species-tree inference

### [ASTRAL](https://pubmed.ncbi.nlm.nih.gov/25161245/)

> Mirarab S, Reaz R, Bayzid MS, Zimmermann T, Swenson MS, Warnow T. ASTRAL: genome-scale coalescent-based species tree estimation. Bioinformatics. 2014;30(17):i541–i548. doi: 10.1093/bioinformatics/btu462. PubMed PMID: 25161245; PubMed Central PMCID: PMC4147915.

> Zhang C, Rabiee M, Sayyari E, Mirarab S. ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees. BMC Bioinformatics. 2018;19(Suppl 6):153. doi: 10.1186/s12859-018-2129-y. PubMed PMID: 29745866; PubMed Central PMCID: PMC6001769.

## Software packaging and containerisation

### [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

> Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018;15(7):475–476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

### [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

> da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017;33(16):2580–2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

### [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

> Merkel D. Docker: lightweight Linux containers for consistent development and deployment. Linux Journal. 2014;2014(239):2. doi: 10.5555/2600239.2600241.

Docker images are used to build or distribute some software environments, although PhyloConstructor executions currently use the Singularity profile.

### [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

> Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. PLoS One. 2017;12(5):e0177459. doi: 10.1371/journal.pone.0177459. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

## Citation recommendations by analysis

For a standard analysis based on downloaded proteomes and concatenation-based phylogenetic inference, cite at least:

* PhyloConstructor
* Nextflow
* the source databases used
* BUSCO
* MAFFT
* trimAl
* PhyKIT
* IQ-TREE
* ModelFinder
* UFBoot2
* ASTRAL or ASTRAL-III
* Singularity
* Bioconda and BioContainers

For an analysis including SRA-derived transcriptomes, additionally cite:

* Sequence Read Archive
* SRA Toolkit
* BBMap
* RNA-SPAdes
* TransDecoder
