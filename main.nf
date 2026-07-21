#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PHYLOCONSTRUCTOR } from './workflows/phyloconstructor'


/*
 * Return a required parameter and stop the pipeline if it is missing
 * or empty.
 */
def requireParameter(String parameterName) {

    def value = params[parameterName]

    if (value == null || value.toString().trim().isEmpty()) {
        error """
Missing required parameter: params.${parameterName}
""".stripIndent()
    }

    return value
}


/*
 * Validate a Boolean parameter.
 *
 * Boolean values must be provided without quotes:
 *
 *     download_uniprot = true
 *
 * and not:
 *
 *     download_uniprot = 'true'
 */
def validateBooleanParameter(String parameterName) {

    def value = params[parameterName]

    if (!(value instanceof Boolean)) {
        error """
Invalid value for params.${parameterName}: ${value}

Expected a Boolean value:

params.${parameterName} = true

or:

params.${parameterName} = false
""".stripIndent()
    }
}


/*
 * Validate the focal taxonomic identifier.
 */
def validateGroupTaxid() {

    def value = requireParameter('group_taxid')
        .toString()
        .trim()

    if (!(value ==~ /^\d+$/)) {
        error """
Invalid value for params.group_taxid: ${value}

The group taxid must contain digits only.
Example:

params.group_taxid = '61985'
""".stripIndent()
    }
}


/*
 * Validate the BUSCO lineage name.
 */
def validateBuscoLineage() {

    def value = requireParameter('busco_lineage')
        .toString()
        .trim()

    if (!(value ==~ /^[A-Za-z0-9._-]+$/)) {
        error """
Invalid value for params.busco_lineage: ${value}

The BUSCO lineage may contain letters, numbers, dots,
underscores and hyphens only.
""".stripIndent()
    }
}


/*
 * Validate the expected number of BUSCO markers.
 */
def validateBuscoDatasetSize() {

    def value = requireParameter('busco_dataset_size')
        .toString()
        .trim()

    if (!(value ==~ /^\d+$/)) {
        error """
Invalid value for params.busco_dataset_size: ${value}

The BUSCO dataset size must be a positive integer.
Example:

params.busco_dataset_size = '1013'
""".stripIndent()
    }

    if ((value as Integer) <= 0) {
        error """
Invalid value for params.busco_dataset_size: ${value}

The BUSCO dataset size must be greater than zero.
""".stripIndent()
    }
}


/*
 * Validate the numeric BUSCO filtering mode provided by the user.
 *
 * Configuration values:
 *
 *   1 = only_single_copy
 *   2 = all_complete
 *   3 = both
 *
 * The function returns the corresponding descriptive internal name.
 */
def validateBuscoFilteringMode() {

    /*
     * Keep the numeric-to-descriptive mapping inside the function.
     *
     * A top-level variable declared with "def" is not directly visible
     * inside script methods in Groovy.
     */
    def buscoFilteringModes = [
        1: 'only_single_copy',
        2: 'all_complete',
        3: 'both'
    ]

    def rawValue = requireParameter('busco_filtering_mode')

    /*
     * Reject Boolean values explicitly because Groovy may otherwise
     * perform unexpected type coercion.
     */
    if (rawValue instanceof Boolean) {
        error """
Invalid value for params.busco_filtering_mode: ${rawValue}

Expected one of the following integer values:

    1 = complete single-copy BUSCO hits only
    2 = all complete BUSCO hits
    3 = test both filtering strategies

Example:

params.busco_filtering_mode = 1
""".stripIndent()
    }

    def modeAsString = rawValue
        .toString()
        .trim()

    if (!(modeAsString ==~ /^[1-3]$/)) {
        error """
Invalid value for params.busco_filtering_mode: ${rawValue}

Expected one of the following integer values:

    1 = complete single-copy BUSCO hits only
    2 = all complete BUSCO hits
    3 = test both filtering strategies

Example:

params.busco_filtering_mode = 1
""".stripIndent()
    }

    def numericMode = modeAsString as Integer

    return buscoFilteringModes[numericMode]
}


/*
 * Return a human-readable description of an internal BUSCO
 * filtering strategy.
 */
def formatBuscoFilteringMode(String mode) {

    switch (mode) {

        case 'only_single_copy':
            return 'complete single-copy BUSCO hits only'

        case 'all_complete':
            return 'all complete BUSCO hits (single-copy and duplicated)'

        case 'both':
            return 'both strategies tested independently'

        default:
            error "Unknown internal BUSCO filtering strategy: ${mode}"
    }
}


/*
 * Read and validate a species CSV file.
 *
 * Expected columns:
 *
 *     taxid,name,fasta
 *
 * The function returns a list of maps with the following keys:
 *
 *     taxid
 *     name
 *     fasta
 */
def readSpeciesCsv(String parameterName) {

    def csvValue = requireParameter(parameterName)
    def csvFile = file(csvValue.toString())

    if (!csvFile.exists()) {
        error """
Input file defined by params.${parameterName} does not exist:

${csvFile}
""".stripIndent()
    }

    if (!csvFile.isFile()) {
        error """
Path defined by params.${parameterName} is not a regular file:

${csvFile}
""".stripIndent()
    }

    if (!csvFile.canRead()) {
        error """
Input file defined by params.${parameterName} is not readable:

${csvFile}
""".stripIndent()
    }

    def lines = csvFile
        .readLines()
        .findAll { line ->
            line != null && !line.trim().isEmpty()
        }

    if (!lines) {
        error """
Input file defined by params.${parameterName} is empty:

${csvFile}
""".stripIndent()
    }

    def header = lines[0]
        .split(',', -1)
        .collect { it.trim().toLowerCase() }

    def expectedHeader = ['taxid', 'name', 'fasta']

    if (header != expectedHeader) {
        error """
Invalid header in params.${parameterName}:

File:
${csvFile}

Observed:
${header.join(',')}

Expected:
${expectedHeader.join(',')}
""".stripIndent()
    }

    def species = []

    lines.drop(1).eachWithIndex { line, index ->

        /*
         * Species names and FASTA paths are not expected to contain commas.
         * Therefore, a simple three-column split is sufficient here.
         */
        def columns = line.split(',', -1)
        def lineNumber = index + 2

        if (columns.size() != 3) {
            error """
Invalid number of columns in params.${parameterName}:

File:
${csvFile}

Line:
${lineNumber}

Content:
${line}

Expected exactly three columns:

taxid,name,fasta
""".stripIndent()
        }

        def taxid = columns[0].trim()
        def name = columns[1].trim()
        def fastaValue = columns[2].trim()

        if (!taxid) {
            error """
Missing taxid in params.${parameterName}:

File:
${csvFile}

Line:
${lineNumber}
""".stripIndent()
        }

        if (!(taxid ==~ /^\d+$/)) {
            error """
Invalid taxid in params.${parameterName}:

File:
${csvFile}

Line:
${lineNumber}

Value:
${taxid}

Taxids must contain digits only.
""".stripIndent()
        }

        if (!name) {
            error """
Missing species name in params.${parameterName}:

File:
${csvFile}

Line:
${lineNumber}
""".stripIndent()
        }

        if (!fastaValue) {
            error """
Missing FASTA path in params.${parameterName}:

File:
${csvFile}

Line:
${lineNumber}

Species:
${name}
""".stripIndent()
        }

        def fastaFile = file(fastaValue)

        if (!fastaFile.exists()) {
            error """
FASTA file listed in params.${parameterName} does not exist:

Species:
${name}

Taxid:
${taxid}

FASTA:
${fastaFile}
""".stripIndent()
        }

        if (!fastaFile.isFile()) {
            error """
FASTA path listed in params.${parameterName} is not a regular file:

Species:
${name}

Taxid:
${taxid}

FASTA:
${fastaFile}
""".stripIndent()
        }

        if (!fastaFile.canRead()) {
            error """
FASTA file listed in params.${parameterName} is not readable:

Species:
${name}

Taxid:
${taxid}

FASTA:
${fastaFile}
""".stripIndent()
        }

        species << [
            taxid: taxid,
            name: name,
            fasta: fastaFile.toAbsolutePath().toString()
        ]
    }

    if (!species) {
        error """
No species were found in params.${parameterName}:

${csvFile}
""".stripIndent()
    }

    /*
     * Taxids must be unique within each input table.
     */
    def duplicatedTaxids = species
        .groupBy { it.taxid }
        .findAll { taxid, records ->
            records.size() > 1
        }

    if (duplicatedTaxids) {
        def details = duplicatedTaxids
            .collect { taxid, records ->
                "${taxid}: ${records.collect { it.name }.join(', ')}"
            }
            .join('\n')

        error """
Duplicated taxids were found in params.${parameterName}:

${details}
""".stripIndent()
    }

    /*
     * Species names must also be unique within each input table.
     */
    def duplicatedNames = species
        .groupBy { it.name.toLowerCase() }
        .findAll { normalizedName, records ->
            records.size() > 1
        }

    if (duplicatedNames) {
        def details = duplicatedNames
            .collect { normalizedName, records ->
                records.collect { it.name }.join(', ')
            }
            .join('\n')

        error """
Duplicated species names were found in params.${parameterName}:

${details}
""".stripIndent()
    }

    return species
}


/*
 * Validate that focal-group and outgroup tables do not contain the
 * same taxids or species names.
 */
def validateSpeciesSets(
    List groupSpecies,
    List outgroupSpecies
) {

    def groupTaxids = groupSpecies.collect { it.taxid }.toSet()
    def outgroupTaxids = outgroupSpecies.collect { it.taxid }.toSet()

    def sharedTaxids = groupTaxids.intersect(outgroupTaxids)

    if (sharedTaxids) {
        def details = sharedTaxids
            .sort()
            .collect { taxid ->

                def groupRecord = groupSpecies.find {
                    it.taxid == taxid
                }

                def outgroupRecord = outgroupSpecies.find {
                    it.taxid == taxid
                }

                "${taxid}: ${groupRecord.name} / ${outgroupRecord.name}"
            }
            .join('\n')

        error """
Taxids are shared between the focal-group and outgroup tables:

${details}
""".stripIndent()
    }

    def groupNames = groupSpecies
        .collect { it.name.toLowerCase() }
        .toSet()

    def outgroupNames = outgroupSpecies
        .collect { it.name.toLowerCase() }
        .toSet()

    def sharedNames = groupNames.intersect(outgroupNames)

    if (sharedNames) {
        def details = sharedNames
            .sort()
            .join('\n')

        error """
Species names are shared between the focal-group and outgroup tables:

${details}
""".stripIndent()
    }
}


/*
 * Validate the combinations of BUSCO completeness and gene occupancy
 * thresholds used for phylogenetic analyses.
 */
def validatePhylogenyThresholds() {

    def thresholds = params.phylogeny_thresholds

    if (thresholds == null) {
        error """
Missing required parameter: params.phylogeny_thresholds

Expected format:

params.phylogeny_thresholds = [
    [60, 80],
    [70, 90],
]
""".stripIndent()
    }

    if (!(thresholds instanceof Collection) || thresholds.isEmpty()) {
        error """
Invalid value for params.phylogeny_thresholds: ${thresholds}

At least one threshold combination must be provided.
Example:

params.phylogeny_thresholds = [
    [60, 80],
]
""".stripIndent()
    }

    def invalidThresholds = thresholds.findAll { threshold ->
        !(threshold instanceof List) ||
        threshold.size() != 2 ||
        !(threshold[0] instanceof Number) ||
        !(threshold[1] instanceof Number)
    }

    if (invalidThresholds) {
        error """
Invalid phylogeny threshold combinations:

${invalidThresholds}

Each combination must contain exactly two numeric values:

[BUSCO completeness, gene occupancy]
""".stripIndent()
    }

    def normalizedThresholds = thresholds.collect { threshold ->
        [
            threshold[0] as Integer,
            threshold[1] as Integer
        ]
    }

    def outOfRangeThresholds = normalizedThresholds.findAll { threshold ->

        def buscoCompleteness = threshold[0]
        def geneOccupancy = threshold[1]

        buscoCompleteness < 0 ||
        buscoCompleteness > 100 ||
        geneOccupancy < 0 ||
        geneOccupancy > 100
    }

    if (outOfRangeThresholds) {
        error """
Phylogeny thresholds must be between 0 and 100:

${outOfRangeThresholds}
""".stripIndent()
    }

    if (normalizedThresholds.size() != normalizedThresholds.unique().size()) {
        error """
Duplicated combinations were found in params.phylogeny_thresholds:

${normalizedThresholds}

Each threshold combination must be unique.
""".stripIndent()
    }
}


/*
 * Validate all pipeline parameters before any process is submitted.
 */
def validateParameters() {

    validateGroupTaxid()

    def groupSpecies = readSpeciesCsv('group_species_csv')
    def outgroupSpecies = readSpeciesCsv('outgroups')

    validateSpeciesSets(
        groupSpecies,
        outgroupSpecies
    )

    [
        'download_uniprot',
        'download_ncbi_proteomes',
        'download_ncbi_assemblies',
        'download_tsa_transcriptomes',
        'download_sra_reads'
    ].each { parameterName ->
        validateBooleanParameter(parameterName)
    }

    validateBuscoLineage()
    validateBuscoDatasetSize()

    /*
     * Convert the user-facing numeric value to the explicit strategy
     * name used internally by the workflow.
     */
    def buscoFilteringStrategy = validateBuscoFilteringMode()

    validatePhylogenyThresholds()

    return [
        groupSpecies          : groupSpecies,
        outgroupSpecies       : outgroupSpecies,
        buscoFilteringStrategy: buscoFilteringStrategy
    ]
}


/*
 * Format a Boolean parameter for the startup summary.
 */
def formatEnabled(Boolean value) {
    value ? 'enabled' : 'disabled'
}


/*
 * Format a species list for the startup summary.
 */
def formatSpeciesList(List species) {

    species
        .sort { a, b ->
            a.name.toLowerCase() <=> b.name.toLowerCase()
        }
        .collect { record ->
            "    - ${record.name} (taxid: ${record.taxid})"
        }
        .join('\n')
}


/*
 * Print the complete validated pipeline configuration before starting
 * any process.
 */
def printParameterSummary(
    List groupSpecies,
    List outgroupSpecies,
    String buscoFilteringStrategy
) {

    def thresholds = params.phylogeny_thresholds
        .collect { threshold ->
            [
                threshold[0] as Integer,
                threshold[1] as Integer
            ]
        }
        .collect { threshold ->
            "    - BUSCO completeness >= ${threshold[0]}%; " +
            "gene occupancy >= ${threshold[1]}%"
        }
        .join('\n')

    def groupSpeciesSummary = formatSpeciesList(groupSpecies)
    def outgroupSpeciesSummary = formatSpeciesList(outgroupSpecies)

    def totalSpecies =
        groupSpecies.size() + outgroupSpecies.size()

    log.info """
===============================================================================
PHYLOCONSTRUCTOR PARAMETERS
===============================================================================

Study group
-----------
Group taxid                 : ${params.group_taxid}
Group species CSV           : ${file(params.group_species_csv).toAbsolutePath()}
Outgroups CSV               : ${file(params.outgroups).toAbsolutePath()}

Species summary
---------------
Focal-group species         : ${groupSpecies.size()}
Outgroup species            : ${outgroupSpecies.size()}
Total user-provided species : ${totalSpecies}

Focal-group species
-------------------
${groupSpeciesSummary}

Outgroup species
----------------
${outgroupSpeciesSummary}

External data retrieval
-----------------------
UniProt proteomes           : ${formatEnabled(params.download_uniprot)}
NCBI proteomes              : ${formatEnabled(params.download_ncbi_proteomes)}
NCBI assemblies             : ${formatEnabled(params.download_ncbi_assemblies)}
TSA transcriptomes          : ${formatEnabled(params.download_tsa_transcriptomes)}
SRA reads                   : ${formatEnabled(params.download_sra_reads)}

BUSCO
-----
Lineage                     : ${params.busco_lineage}
Expected dataset size       : ${params.busco_dataset_size}
Filtering mode              : ${params.busco_filtering_mode}
Internal filtering strategy : ${buscoFilteringStrategy}
Filtering description       : ${formatBuscoFilteringMode(
                                    buscoFilteringStrategy
                                )}

Phylogenetic analyses
---------------------
Threshold combinations      :
${thresholds}

Execution
---------
Executor                    : SLURM
Nextflow profile            : ${workflow.profile ?: 'default'}
Nextflow version            : ${workflow.nextflow.version}
Project directory           : ${projectDir}
Launch directory            : ${launchDir}
Work directory              : ${workDir}

===============================================================================
""".stripIndent()
}


workflow {

    /*
     * Validate the complete configuration before submitting processes.
     */
    def validatedInputs = validateParameters()

    /*
     * Print the original numeric mode and the corresponding internal
     * BUSCO filtering strategy.
     */
    printParameterSummary(
        validatedInputs.groupSpecies,
        validatedInputs.outgroupSpecies,
        validatedInputs.buscoFilteringStrategy
    )

    /*
     * From this point onward, workflow components only receive the
     * descriptive internal strategy name.
     *
     * Examples:
     *   only_single_copy
     *   all_complete
     *   both
     */

    PHYLOCONSTRUCTOR(
        validatedInputs.buscoFilteringStrategy
    )
}