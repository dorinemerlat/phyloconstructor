#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})


# ============================================================
# Input / output
# ============================================================

TREE_FILE <- "/shared/projects/metainvert/phyloconstructor2/cache/phylogeny/render_tree_iqtree/busco-c_60_gene-occupancy_80/trimal_only_single_copy_busco-c_60_gene-occupancy_80_iqtree_taxonomic_tree_without_clade_annotations.nwk"

CALIBRATION_FILE <- "/shared/projects/metainvert/phyloconstructor2/calibration_benavides.csv"

OUTPUT_TREE <- "/shared/projects/metainvert/phyloconstructor2/cache/cafe/tree_renamed_iqtree_pruned_dated_all.nwk"

OUTPUT_CALIBRATIONS <- "/shared/projects/metainvert/phyloconstructor2/cache/cafe/chronos_calibrations_used_all.csv"

OUTPUT_NODE_AGES <- "/shared/projects/metainvert/phyloconstructor2/cache/cafe/chronos_node_ages_all.csv"

OUTPUT_PLOT <- "/shared/projects/metainvert/phyloconstructor2/cache/cafe/tree_renamed_iqtree_pruned_dated_all.pdf"


# ============================================================
# chronos parameters
# ============================================================

CHRONOS_LAMBDA <- 1
CHRONOS_MODEL <- "correlated"


# ============================================================
# Read input tree
# ============================================================

cat("Reading tree:\n")
cat("  ", TREE_FILE, "\n\n")


# Read raw Newick text.
newick <- paste(
  readLines(TREE_FILE, warn = FALSE),
  collapse = ""
)

# Remove empty internal node labels produced by ETE3:
# )"" -> )
newick <- gsub('\\)""', ')', newick)

# Convert double-quoted taxon names to standard Newick single quotes:
# "Argiope bruennichi" -> 'Argiope bruennichi'
newick <- gsub(
  '"([^"]+)"',
  "'\\1'",
  newick
)

# Read the cleaned Newick.
tree <- read.tree(text = newick)
tree$tip.label <- gsub("^'|'$", "", tree$tip.label)

cat(
  "Argiope found:",
  "Argiope bruennichi" %in% tree$tip.label,
  "\n"
)
if (is.null(tree)) {
  stop("Could not read input tree.")
}

cat("Tips :", Ntip(tree), "\n")
cat("Nodes:", Nnode(tree), "\n")

if (!is.rooted(tree)) {
  stop(
    "The input tree is not rooted. ",
    "Root the IQ-TREE tree before running chronos."
  )
}

if (is.null(tree$edge.length)) {
  stop("The input tree has no branch lengths.")
}

if (any(tree$edge.length <= 0)) {
  stop(
    "The input tree contains zero or negative branch lengths. ",
    "chronos requires positive molecular branch lengths."
  )
}

# IQ-TREE bootstrap/support labels are not needed for dating.
tree$node.label <- NULL

# Remove an optional root edge.
tree$root.edge <- NULL


# ============================================================
# Read Benavides et al. calibration table
# ============================================================

cat("\nReading calibration table:\n")
cat("  ", CALIBRATION_FILE, "\n\n")

benavides <- read.csv(
  CALIBRATION_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "id",
  "crown_group_node",
  "minimum_age",
  "mcmctree_bound"
)

missing_columns <- setdiff(
  required_columns,
  colnames(benavides)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns in calibration CSV: ",
    paste(missing_columns, collapse = ", ")
  )
}

benavides$id <- as.integer(benavides$id)


# ============================================================
# Helper functions
# ============================================================

parse_minimum_age <- function(x) {
  if (is.na(x) || x == "") {
    return(NA_real_)
  }
  
  value <- sub(
    "^\\s*([0-9.]+).*$",
    "\\1",
    x
  )
  
  as.numeric(value)
}


parse_mcmctree_bound <- function(x) {
  if (is.na(x) || x == "") {
    return(
      list(
        type = NA_character_,
        min = NA_real_,
        max = NA_real_
      )
    )
  }
  
  x <- gsub("\\s+", "", x)
  
  # B(4.76,5.28)
  if (grepl("^B\\(", x)) {
    values <- sub(
      "^B\\(([^,]+),([^\\)]+)\\)$",
      "\\1,\\2",
      x
    )
    
    values <- as.numeric(
      strsplit(values, ",")[[1]]
    )
    
    return(
      list(
        type = "B",
        min = values[1] * 100,
        max = values[2] * 100
      )
    )
  }
  
  # U(5.42)
  if (grepl("^U\\(", x)) {
    value <- sub(
      "^U\\(([^\\)]+)\\)$",
      "\\1",
      x
    )
    
    return(
      list(
        type = "U",
        min = 0,
        max = as.numeric(value) * 100
      )
    )
  }
  
  stop(
    "Unsupported MCMCTree bound format: ",
    x
  )
}


get_calibration_ages <- function(calibration_id) {
  row <- benavides[
    benavides$id == calibration_id,
    ,
    drop = FALSE
  ]
  
  if (nrow(row) != 1) {
    stop(
      "Expected exactly one row for calibration ID ",
      calibration_id,
      ", found ",
      nrow(row)
    )
  }
  
  # Special adaptation for the root.
  #
  # Benavides et al. list Rusophycus at 528 Ma for the root,
  # while the base of the Cambrian (542 Ma) is used as a deep
  # maximum elsewhere in their calibration scheme.
  #
  # For chronos, we use the interval 528–542 Ma.
  if (calibration_id == 0) {
    return(
      list(
        min = 528,
        max = 542
      )
    )
  }
  
  fossil_min <- parse_minimum_age(
    row$minimum_age
  )
  
  bound <- parse_mcmctree_bound(
    row$mcmctree_bound
  )
  
  # Preserve the more precise fossil minimum from the table
  # rather than the rounded lower value in B(...).
  age_min <- if (
    !is.na(fossil_min)
  ) {
    fossil_min
  } else {
    bound$min
  }
  
  list(
    min = age_min,
    max = bound$max
  )
}


check_tip <- function(tip) {
  if (!tip %in% tree$tip.label) {
    stop(
      "Species not found in tree: ",
      tip
    )
  }
}


get_mrca_node <- function(tip1, tip2) {
  check_tip(tip1)
  check_tip(tip2)
  
  node <- getMRCA(
    tree,
    c(tip1, tip2)
  )
  
  if (is.null(node) || is.na(node)) {
    stop(
      "Could not identify MRCA for: ",
      tip1,
      " / ",
      tip2
    )
  }
  
  node
}


# ============================================================
# Calibrations applicable to THIS pruned tree
# ============================================================
#
# The IDs correspond to Table 2 of Benavides et al. (2023).
#
# We deliberately do not use calibrations for clades whose
# defining lineages are absent from this pruned tree.
#
# Examples:
#   Crown Chilopoda:
#       requires Scutigeromorpha + Pleurostigmophora
#       -> Scutigeromorpha absent
#
#   Crown Diplopoda:
#       requires Penicillata + Chilognatha
#       -> Penicillata absent
#
# ============================================================

# ============================================================
# Calibrations applicable to the complete 141-species tree
# ============================================================
#
# IDs correspond to Table 2 of Benavides et al. (2023).
#
# The two tips are representatives of the two descendant
# lineages defining the calibrated crown node.
#
# Some species used by Benavides et al. are absent from this
# tree. In those cases, another species from the same relevant
# descendant lineage is used.
#
# Calibrations are excluded when the required lineage is absent
# and the corresponding crown node therefore cannot be located
# unambiguously.
# ============================================================

calibration_map <- data.frame(
  id = c(
    0,   # Root
    2,   # Crown Euchelicerata (partim)
    3,   # Crown Mandibulata
    4,   # Crown Pancrustacea
    6,   # Crown Allotriocarida (partim)
    7,   # Crown Pauropoda
    8,   # Crown Symphyla
    10,  # Crown Chilopoda
    12,  # Crown Pleurostigmophora
    13,  # Crown Lithobiomorpha
    14,  # Crown Epimorpha
    16,  # Crown Scolopendromorpha
    17,  # Crown Tykhepoda
    18,  # Crown Scolopocryptopidae
    19,  # Crown Diplopoda
    20,  # Crown Polyxenidae
    21,  # Crown Chilognatha
    23,  # Crown Helminthomorpha
    24,  # Crown Platydesmida
    25,  # Crown Polyzoniida
    27,  # Crown Coelochaeta
    29,  # Crown Polydesmida (partim)
    30,  # Crown Spirostreptida
    31   # Crown Blaniuidea + Nemasomatoidea + Juloidea
  ),
  
  tip1 = c(
    NA,
    
    # 2 Crown Euchelicerata
    "Tachypleus tridentatus",
    
    # 3 Crown Mandibulata
    "Daphnia pulex",
    
    # 4 Crown Pancrustacea
    "Eriocheir sinensis",
    
    # 6 Crown Allotriocarida
    "Daphnia pulex",
    
    # 7 Crown Pauropoda
    "Pauropus huxleyi",
    
    # 8 Crown Symphyla
    "Symphylella vulgaris",
    
    # 10 Crown Chilopoda
    "Scutigera coleoptrata",
    
    # 12 Crown Pleurostigmophora
    "Craterostigmus tasmanianus",
    
    # 13 Crown Lithobiomorpha
    "Anopsobius giribeti",
    
    # 14 Crown Epimorpha
    "Strigamia maritima",
    
    # 16 Crown Scolopendromorpha
    "Scolopendra cretica",
    
    # 17 Crown Tykhepoda
    "Cryptops parisi",
    
    # 18 Crown Scolopocryptopidae
    "Newportia adisi",
    
    # 19 Crown Diplopoda
    "Polyxenus lagurus",
    
    # 20 Crown Polyxenidae
    "Polyxenus lagurus",
    
    # 21 Crown Chilognatha
    "Glomeris marginata",
    
    # 23 Crown Helminthomorpha
    "Brachycybe producta",
    
    # 24 Crown Platydesmida
    "Andrognathus corticarius",
    
    # 25 Crown Polyzoniida
    "Rhinotus purpureus",
    
    # 27 Crown Coelochaeta
    "Abacion magnum",
    
    # 29 Crown Polydesmida partim
    "Polydesmus complanatus",
    
    # 30 Crown Spirostreptida
    "Cambala annulata",
    
    # 31 Crown Blaniuidea + Nemasomatoidea + Juloidea
    "Uroblaniulus sp. MITS883"
  ),
  
  tip2 = c(
    NA,
    
    # 2 Crown Euchelicerata
    "Argiope bruennichi",
    
    # 3 Crown Mandibulata
    "Scutigera coleoptrata",
    
    # 4 Crown Pancrustacea
    "Daphnia pulex",
    
    # 5 Crown Multicrustacea
    "Eriocheir sinensis",
    
    # 6 Crown Allotriocarida
    "Folsomia candida",
    
    # 7 Crown Pauropoda
    "Acopauropus ornatus",
    
    # 8 Crown Symphyla
    "Hanseniella nivea",
    
    # 10 Crown Chilopoda
    "Lithobius variegatus",
    
    # 12 Crown Pleurostigmophora
    "Lithobius variegatus",
    
    # 13 Crown Lithobiomorpha
    "Lithobius variegatus",
    
    # 14 Crown Epimorpha
    "Scolopendra cretica",
    
    # 16 Crown Scolopendromorpha
    "Cryptops parisi",
    
    # 17 Crown Tykhepoda
    "Theatops posticus",
    
    # 18 Crown Scolopocryptopidae
    "Scolopocryptops rubiginosus",
    
    # 19 Crown Diplopoda
    "Glomeris marginata",
    
    # 20 Crown Polyxenidae
    "Eudigraphis taiwaniensis",
    
    # 21 Crown Chilognatha
    "Brachycybe producta",
    
    # 23 Crown Helminthomorpha
    "Polydesmus complanatus",
    
    # 24 Crown Platydesmida
    "Platydesmus sp. MITS438",
    
    # 25 Crown Polyzoniida
    "Polyzonium germanicum",
    
    # 27 Crown Coelochaeta
    "Chordeuma sylvestre",
    
    # 29 Crown Polydesmida partim
    "Cyrtodesmus sp. MITS521",
    
    # 30 Crown Spirostreptida
    "Archispirostreptus syriacus",
    
    # 31 Crown Blaniuidea + Nemasomatoidea + Juloidea
    "Julus scandinavius"
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# Identify calibrated nodes
# ============================================================

root_node <- Ntip(tree) + 1

calibration_results <- list()

for (i in seq_len(nrow(calibration_map))) {
  
  calibration_id <- calibration_map$id[i]
  
  source_row <- benavides[
    benavides$id == calibration_id,
    ,
    drop = FALSE
  ]
  
  if (nrow(source_row) != 1) {
    stop(
      "Calibration ID ",
      calibration_id,
      " not found uniquely in CSV."
    )
  }
  
  if (calibration_id == 0) {
    node <- root_node
    
    tip1 <- NA_character_
    tip2 <- NA_character_
  } else {
    tip1 <- calibration_map$tip1[i]
    tip2 <- calibration_map$tip2[i]
    
    node <- get_mrca_node(
      tip1,
      tip2
    )
  }
  
  ages <- get_calibration_ages(
    calibration_id
  )
  
  calibration_results[[i]] <- data.frame(
    id = calibration_id,
    crown_group_node = source_row$crown_group_node,
    calibration_fossil = source_row$calibration_fossil,
    tip1 = tip1,
    tip2 = tip2,
    node = node,
    age.min = ages$min,
    age.max = ages$max,
    stringsAsFactors = FALSE
  )
}

calibration_table <- do.call(
  rbind,
  calibration_results
)


# ============================================================
# Check for duplicated calibrated nodes
# ============================================================

duplicated_nodes <- unique(
  calibration_table$node[
    duplicated(calibration_table$node)
  ]
)

if (length(duplicated_nodes) > 0) {
  cat(
    "\nWARNING: multiple calibrations map to the same node:\n"
  )
  
  for (node in duplicated_nodes) {
    print(
      calibration_table[
        calibration_table$node == node,
        ,
        drop = FALSE
      ]
    )
  }
  
  stop(
    "Duplicated calibrated nodes detected. ",
    "Review the calibration mapping before running chronos."
  )
}


# ============================================================
# Display selected calibrations
# ============================================================

cat("\nCalibrations used:\n\n")

print(
  calibration_table[
    ,
    c(
      "id",
      "crown_group_node",
      "node",
      "age.min",
      "age.max",
      "tip1",
      "tip2"
    )
  ],
  row.names = FALSE
)


# ============================================================
# Save calibration mapping
# ============================================================

write.csv(
  calibration_table,
  OUTPUT_CALIBRATIONS,
  row.names = FALSE,
  quote = TRUE
)

cat(
  "\nCalibration table written to:\n  ",
  OUTPUT_CALIBRATIONS,
  "\n"
)


# ============================================================
# Build chronos calibration data frame
# ============================================================

chronos_calibration <- makeChronosCalib(
  tree,
  node = calibration_table$node,
  age.min = calibration_table$age.min,
  age.max = calibration_table$age.max,
  soft.bounds = FALSE
)

cat("\nchronos calibration data frame:\n\n")

print(
  chronos_calibration,
  row.names = FALSE
)


# ============================================================
# Run chronos
# ============================================================

cat("\nRunning chronos...\n")
cat("Model :", CHRONOS_MODEL, "\n")
cat("Lambda:", CHRONOS_LAMBDA, "\n\n")

dated_tree <- chronos(
  tree,
  lambda = CHRONOS_LAMBDA,
  model = CHRONOS_MODEL,
  calibration = chronos_calibration,
  control = chronos.control(
    iter.max = 100000,
    eval.max = 100000
  ),
  quiet = FALSE
)


# ============================================================
# Check ultrametricity
# ============================================================

ultrametric <- is.ultrametric(
  dated_tree
)

cat(
  "\nUltrametric tree: ",
  ultrametric,
  "\n",
  sep = ""
)

root_to_tip <- node.depth.edgelength(
  dated_tree
)[seq_len(Ntip(dated_tree))]

cat(
  "Root-to-tip minimum: ",
  min(root_to_tip),
  " Ma\n",
  sep = ""
)

cat(
  "Root-to-tip maximum: ",
  max(root_to_tip),
  " Ma\n",
  sep = ""
)

cat(
  "Root-to-tip difference: ",
  max(root_to_tip) - min(root_to_tip),
  " Ma\n",
  sep = ""
)

if (!ultrametric) {
  stop(
    "chronos returned a tree that is not ultrametric."
  )
}


# ============================================================
# Calculate node ages
# ============================================================

depths <- node.depth.edgelength(
  dated_tree
)

root_age <- max(
  depths[seq_len(Ntip(dated_tree))]
)

node_numbers <- seq_len(
  Ntip(dated_tree) + Nnode(dated_tree)
)

node_ages <- root_age - depths[node_numbers]

node_type <- ifelse(
  node_numbers <= Ntip(dated_tree),
  "tip",
  "internal"
)

node_name <- rep(
  "",
  length(node_numbers)
)

node_name[
  seq_len(Ntip(dated_tree))
] <- dated_tree$tip.label

node_age_table <- data.frame(
  node = node_numbers,
  type = node_type,
  name = node_name,
  age_Ma = node_ages,
  stringsAsFactors = FALSE
)

write.csv(
  node_age_table,
  OUTPUT_NODE_AGES,
  row.names = FALSE,
  quote = TRUE
)

cat(
  "\nNode ages written to:\n  ",
  OUTPUT_NODE_AGES,
  "\n"
)


# ============================================================
# Write dated tree
# ============================================================

write.tree(
  dated_tree,
  file = OUTPUT_TREE
)

cat(
  "\nDated ultrametric tree written to:\n  ",
  OUTPUT_TREE,
  "\n"
)


# ============================================================
# Plot dated tree
# ============================================================

pdf(
  OUTPUT_PLOT,
  width =6,
  height = 15
)

plot(
  dated_tree,
  type = "phylogram",
  show.tip.label = TRUE,
  cex = 0.45,
  no.margin = FALSE
)

title(
  main = paste0(
    "Time-calibrated tree — chronos\n",
    "Root age = ",
    round(root_age, 2),
    " Ma"
  )
)

dev.off()

cat(
  "Dated tree plot written to:\n  ",
  OUTPUT_PLOT,
  "\n"
)


# ============================================================
# Final summary
# ============================================================

cat("\n========================================\n")
cat("chronos dating completed\n")
cat("========================================\n")
cat("Input tree       :", TREE_FILE, "\n")
cat("Calibration CSV  :", CALIBRATION_FILE, "\n")
cat("Calibrations used:", nrow(calibration_table), "\n")
cat("Model            :", CHRONOS_MODEL, "\n")
cat("Lambda           :", CHRONOS_LAMBDA, "\n")
cat("Root age         :", round(root_age, 6), "Ma\n")
cat("Ultrametric      :", ultrametric, "\n")
cat("Output tree      :", OUTPUT_TREE, "\n")
cat("========================================\n")
