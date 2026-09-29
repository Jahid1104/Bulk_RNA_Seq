# ============================================================
# TOMATO GO ANALYSIS
# ITAG4.0 GO TERMS + GO OBO + MULTIPLE DESEQ2 FILES
# ============================================================
#
# OUTPUT:
#
# 1. Individual GO enrichment results
# 2. Combined GO enrichment results
# 3. Top 30 GO enrichment terms
# 4. Functional annotation results
# 5. Functional annotation figure
# 6. Functional enrichment figure
# 7. CSV files containing plotting data
# 8. 1000 dpi TIFF figures
#
# FIGURE 1:
# Functional annotation
#   - 8 rows = comparisons (a-h, letters on the right)
#   - 3 columns = Molecular Function, Cellular Component, Biological Process
#   - Vertical bars, GO term labels read vertically (parallel to the bars)
#   - UP = red
#   - DOWN = steelblue
#   - Numbers = UP/DOWN genes
#
# FIGURE 2:
# Functional enrichment
#   - Top 30 GO terms (top 10 per category, grouped MF / CC / BP)
#   - 8 comparison columns (a-h)
#   - Horizontal bars
#   - Length = number of genes (UP + DOWN combined)
#   - Color = -log10(FDR/Q value)
#   - Number = gene count
#
# All input/output paths are hardcoded below (no file-choose dialogs).
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "purrr",
  "tibble",
  "ggtext"
)

for (p in cran_packages) {
  
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  
}


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


if (!requireNamespace("topGO", quietly = TRUE)) {
  
  BiocManager::install(
    "topGO",
    ask = FALSE,
    update = FALSE
  )
  
}


library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(topGO)
library(ggtext)


# ------------------------------------------------------------
# topGO (via AnnotationDbi/Biobase) defines its own select()
# and filter(), which mask dplyr's versions after loading.
# Force dplyr's versions back so the rest of the script works
# as written.
# ------------------------------------------------------------

select <- dplyr::select

filter <- dplyr::filter


# ============================================================
# 2. USER SETTINGS
# ============================================================

# ------------------------------------------------------------
# DEG thresholds
# ------------------------------------------------------------

PADJ_CUTOFF <- 0.05

LFC_CUTOFF <- 1


# ------------------------------------------------------------
# Number of GO terms for functional annotation
# ------------------------------------------------------------

TOP_ANNOTATION_TERMS <- 10


# ------------------------------------------------------------
# Number of GO terms for functional enrichment
# 10 per category (MF, CC, BP) = 30 terms total
# ------------------------------------------------------------

TOP_ENRICHMENT_TERMS_PER_CATEGORY <- 10

TOP_ENRICHMENT_TERMS <- TOP_ENRICHMENT_TERMS_PER_CATEGORY * 3


# ------------------------------------------------------------
# Publication figure font sizes
# ------------------------------------------------------------
#
# Increase these if needed.
#
# ------------------------------------------------------------

FONT_BASE <- 13

FONT_AXIS_X <- 12

FONT_AXIS_Y <- 13

FONT_STRIP <- 16

FONT_LEGEND <- 14

FONT_GO <- 12   # increased from 9 -- GO term/name labels on the functional annotation figure (now bold)

FONT_NUMBER <- 5   # UP/DOWN count labels on the functional annotation figure (now printed vertically, parallel to the GO term labels)

FONT_TITLE <- 14


# ------------------------------------------------------------
# Functional enrichment figure fonts/sizing
# Matched directly to the KEGG functional enrichment figure
# (6_2_KEGG_Analysis.R, Figure 2) so the two look like a pair.
# ------------------------------------------------------------

FONT_BASE_ENRICH    <- 9    # overall base size
FONT_AXIS_ENRICH    <- 9    # numeric x-axis (gene count) text
FONT_GO_ENRICH      <- 10   # GO term labels (bold), same role as KEGG's FONT_PATHWAY
FONT_STRIP_ENRICH   <- 12   # comparison a-h strip labels
FONT_LEGEND_ENRICH  <- 9
FONT_NUMBER_ENRICH  <- 2.6  # gene-count data labels at the end of each bar


# ------------------------------------------------------------
# Figure dimensions
# ------------------------------------------------------------

ANNOTATION_WIDTH <- 18

ANNOTATION_HEIGHT <- 22

ENRICHMENT_WIDTH <- 18   # matches KEGG Figure 2 (FIG2_WIDTH)

ENRICHMENT_HEIGHT <- 12   # matches KEGG Figure 2 (FIG2_HEIGHT)


# ------------------------------------------------------------
# TIFF resolution
# ------------------------------------------------------------

FIG_DPI <- 1000


# ------------------------------------------------------------
# GO category display order (left to right / top to bottom)
# ------------------------------------------------------------
#
# Molecular Function, Cellular Component, Biological Process
#
# ------------------------------------------------------------

ONTOLOGY_LEVELS <- c("MF", "CC", "BP")

ONTOLOGY_LABELS <- c(
  MF = "Molecular Function",
  CC = "Cellular Component",
  BP = "Biological Process"
)


# ------------------------------------------------------------
# Shorten "biological process" / "cellular component" /
# "molecular function" to BP / CC / MF wherever they appear
# inside a GO term name, to keep labels short.
# ------------------------------------------------------------

shorten_go_name <- function(x) {
  
  x <- str_replace_all(
    x,
    regex("biological_process|biological process", ignore_case = TRUE),
    "BP"
  )
  
  x <- str_replace_all(
    x,
    regex("cellular_component|cellular component", ignore_case = TRUE),
    "CC"
  )
  
  x <- str_replace_all(
    x,
    regex("molecular_function|molecular function", ignore_case = TRUE),
    "MF"
  )
  
  x
  
}


# ------------------------------------------------------------
# Trim a GO term name down to its first 3 words, adding "..."
# if anything was cut off, to keep labels short.
# ------------------------------------------------------------

truncate_go_name <- function(x, max_words = 3) {
  
  words <- strsplit(x, "\\s+")
  
  vapply(
    words,
    function(w) {
      
      if (length(w) > max_words) {
        
        paste0(
          paste(w[seq_len(max_words)], collapse = " "),
          " ..."
        )
        
      } else {
        
        paste(w, collapse = " ")
        
      }
      
    },
    character(1)
  )
  
}


# ------------------------------------------------------------
# Colors used to tell MF / CC / BP terms apart on the
# functional enrichment figure.
# ------------------------------------------------------------

ONTOLOGY_COLORS <- c(
  MF = "#1b9e77",
  CC = "#d95f02",
  BP = "#7570b3"
)


# ------------------------------------------------------------
# Build the two-line GO term label (name + GO ID) used on
# both figures. `colored = TRUE` wraps the label in a markdown
# color span (functional enrichment figure only), rendered via
# ggtext::element_markdown().
# ------------------------------------------------------------

make_go_label <- function(GO_name, GO_ID, Ontology, colored = FALSE) {
  
  Ontology <- as.character(Ontology)
  
  name_txt <-
    truncate_go_name(
      shorten_go_name(GO_name)
    )
  
  if (colored) {
    
    col <- ONTOLOGY_COLORS[Ontology]
    
    paste0(
      "<span style='color:", col, "'>",
      name_txt,
      "<br>",
      GO_ID,
      "</span>"
    )
    
  } else {
    
    paste0(
      name_txt,
      "\n",
      GO_ID
    )
    
  }
  
}


# ============================================================
# 3. ITAG GO FILE, OBO FILE, DESEQ2 FILES, OUTPUT FOLDER
# ============================================================
#
# All paths are fixed below. Edit these if your file locations
# change.
#
# ============================================================

go_file <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/Reference/ITAG4.0_goterms.txt"

obo_file <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/Reference/go-basic.obo"

deseq_dir <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/New_Name/Results/5_DESeq2/Tables/Comparisons"

output_dir <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/New_Name/Results/6_GO"


# ------------------------------------------------------------
# DESeq2 comparison files, in letter order a - h
# ------------------------------------------------------------

deseq_filenames <- c(
  a = "CLN1466EA vs NC123S_DESeq2_results.csv",
  b = "Flower vs Leaf_DESeq2_results.csv",
  c = "CLN1466EA_Flower vs NC123S_Flower_DESeq2_results.csv",
  d = "CLN1466EA_Flower_24h vs NC123S_Flower_24h_DESeq2_results.csv",
  e = "CLN1466EA_Flower_72h vs NC123S_Flower_72h_DESeq2_results.csv",
  f = "CLN1466EA_Leaf vs NC123S_Leaf_DESeq2_results.csv",
  g = "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h_DESeq2_results.csv",
  h = "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h_DESeq2_results.csv"
)

deseq_files <- file.path(deseq_dir, deseq_filenames)

names(deseq_files) <- names(deseq_filenames)


if (!file.exists(go_file)) {
  stop("ITAG GO terms file not found:\n", go_file)
}

if (!file.exists(obo_file)) {
  stop("OBO file not found:\n", obo_file)
}

missing_deseq <- deseq_files[!file.exists(deseq_files)]

if (length(missing_deseq) > 0) {
  stop(
    "The following DESeq2 file(s) were not found:\n",
    paste(missing_deseq, collapse = "\n")
  )
}


cat(
  "\nNumber of DESeq2 files:",
  length(deseq_files),
  "\n"
)


# ============================================================
# 4. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


individual_dir <-
  
  file.path(
    output_dir,
    "Individual_Results"
  )


dir.create(
  individual_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


figure_dir <-
  
  file.path(
    output_dir,
    "Figures"
  )


dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 5. READ ITAG GO FILE
# ============================================================
#
# Your file has lines like:
#
# Solyc01g005000.2<TAB>GO:0016831,GO:0019752,GO:0030170
#
# Some genes have no GO term and therefore contain only:
#
# Solyc01g004000.1
#
# We read the file line-by-line.
#
# ============================================================

cat("\nReading ITAG GO file...\n")


raw_lines <-
  
  readLines(
    go_file,
    warn = FALSE
  )


raw_lines <-
  
  raw_lines[
    nzchar(
      trimws(raw_lines)
    )
  ]


# Remove comment lines

raw_lines <-
  
  raw_lines[
    !grepl(
      "^#",
      raw_lines
    )
  ]


gene_vector <-
  character(
    length(raw_lines)
  )


go_vector <-
  character(
    length(raw_lines)
  )


for (i in seq_along(raw_lines)) {
  
  pieces <-
    
    strsplit(
      raw_lines[i],
      "\t",
      fixed = TRUE
    )[[1]]
  
  
  gene_vector[i] <-
    
    trimws(
      pieces[1]
    )
  
  
  if (length(pieces) >= 2) {
    
    go_vector[i] <-
      
      trimws(
        paste(
          pieces[-1],
          collapse = "\t"
        )
      )
    
  } else {
    
    go_vector[i] <- ""
    
  }
  
}


go_raw <-
  
  tibble(
    
    Gene =
      gene_vector,
    
    GO =
      go_vector
    
  )


# ============================================================
# 6. CLEAN GENE IDs
# ============================================================

go_raw$Gene <-
  
  str_extract(
    go_raw$Gene,
    "Solyc[0-9]{2}g[0-9]+\\.[0-9]+"
  )


go_raw <-
  
  go_raw %>%
  
  filter(
    !is.na(Gene)
  )


# ============================================================
# 7. CREATE gene2GO
# ============================================================

gene2GO <- list()


for (i in seq_len(nrow(go_raw))) {
  
  gene <- go_raw$Gene[i]
  
  go_string <- go_raw$GO[i]
  
  
  if (
    is.na(go_string) ||
    go_string == ""
  ) {
    next
  }
  
  
  terms <-
    
    unlist(
      strsplit(
        go_string,
        ","
      )
    )
  
  
  terms <-
    
    trimws(
      terms
    )
  
  
  terms <-
    
    terms[
      grepl(
        "^GO:[0-9]{7}$",
        terms
      )
    ]
  
  
  if (length(terms) == 0) {
    next
  }
  
  
  gene2GO[[gene]] <-
    
    unique(
      terms
    )
  
}


gene2GO <-
  
  gene2GO[
    lengths(gene2GO) > 0
  ]


cat(
  "Genes with GO annotations:",
  length(gene2GO),
  "\n"
)


if (length(gene2GO) == 0) {
  
  stop(
    "No valid gene-GO associations were found."
  )
  
}


# ============================================================
# 8. READ OBO FILE
# ============================================================

cat(
  "\nReading OBO file...\n"
)


obo_lines <-
  
  readLines(
    obo_file,
    warn = FALSE
  )


term_start <-
  
  which(
    trimws(
      obo_lines
    ) == "[Term]"
  )


obo_list <- list()


for (i in seq_along(term_start)) {
  
  start <- term_start[i]
  
  
  if (i < length(term_start)) {
    
    end <-
      term_start[i + 1] - 1
    
  } else {
    
    end <-
      length(obo_lines)
    
  }
  
  
  block <-
    obo_lines[start:end]
  
  
  id_line <-
    
    grep(
      "^id:",
      block,
      value = TRUE
    )
  
  
  name_line <-
    
    grep(
      "^name:",
      block,
      value = TRUE
    )
  
  
  namespace_line <-
    
    grep(
      "^namespace:",
      block,
      value = TRUE
    )
  
  
  if (
    length(id_line) == 0 ||
    length(name_line) == 0 ||
    length(namespace_line) == 0
  ) {
    next
  }
  
  
  go_id <-
    
    sub(
      "^id:\\s*",
      "",
      id_line[1]
    )
  
  
  go_name <-
    
    sub(
      "^name:\\s*",
      "",
      name_line[1]
    )
  
  
  namespace <-
    
    sub(
      "^namespace:\\s*",
      "",
      namespace_line[1]
    )
  
  
  # Remove obsolete GO terms
  
  if (
    any(
      grepl(
        "^is_obsolete:\\s*true",
        block
      )
    )
  ) {
    next
  }
  
  
  obo_list[[
    length(obo_list) + 1
  ]] <-
    
    tibble(
      
      GO_ID =
        go_id,
      
      GO_name =
        go_name,
      
      namespace =
        namespace
      
    )
  
}


go_info <-
  
  bind_rows(
    obo_list
  )


# ============================================================
# 9. BROAD GO CATEGORIES
# ============================================================

go_info <-
  
  go_info %>%
  
  mutate(
    
    Ontology = case_when(
      
      namespace ==
        "biological_process" ~
        "BP",
      
      namespace ==
        "cellular_component" ~
        "CC",
      
      namespace ==
        "molecular_function" ~
        "MF",
      
      TRUE ~
        NA_character_
      
    )
    
  ) %>%
  
  filter(
    !is.na(Ontology)
  )


cat(
  "GO terms:",
  nrow(go_info),
  "\n"
)


# ============================================================
# 10. SPLIT gene2GO BY ONTOLOGY
# ============================================================

gene2GO_BP <-
  
  lapply(
    gene2GO,
    function(x) {
      
      intersect(
        x,
        go_info$GO_ID[
          go_info$Ontology == "BP"
        ]
      )
      
    }
  )


gene2GO_CC <-
  
  lapply(
    gene2GO,
    function(x) {
      
      intersect(
        x,
        go_info$GO_ID[
          go_info$Ontology == "CC"
        ]
      )
      
    }
  )


gene2GO_MF <-
  
  lapply(
    gene2GO,
    function(x) {
      
      intersect(
        x,
        go_info$GO_ID[
          go_info$Ontology == "MF"
        ]
      )
      
    }
  )


gene2GO_BP <-
  gene2GO_BP[
    lengths(gene2GO_BP) > 0
  ]


gene2GO_CC <-
  gene2GO_CC[
    lengths(gene2GO_CC) > 0
  ]


gene2GO_MF <-
  gene2GO_MF[
    lengths(gene2GO_MF) > 0
  ]


# ============================================================
# 11. READ DESEQ2 RESULT
# ============================================================

read_deseq <- function(file) {
  
  
  cat(
    "\nReading:",
    basename(file),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Try tab-delimited first
  # ----------------------------------------------------------
  
  dat <-
    
    tryCatch(
      
      read.delim(
        file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      
      error =
        function(e) NULL
      
    )
  
  
  # ----------------------------------------------------------
  # If only one column, try comma
  # ----------------------------------------------------------
  
  if (
    is.null(dat) ||
    ncol(dat) <= 1
  ) {
    
    dat <-
      
      tryCatch(
        
        read.csv(
          file,
          check.names = FALSE,
          stringsAsFactors = FALSE
        ),
        
        error =
          function(e) NULL
        
      )
    
  }
  
  
  if (
    is.null(dat) ||
    ncol(dat) <= 1
  ) {
    
    stop(
      paste0(
        "\nCould not properly read:\n",
        file,
        "\n\nCheck whether the file is tab- or comma-separated."
      )
    )
    
  }
  
  
  cat(
    "Columns:",
    paste(
      colnames(dat),
      collapse = ", "
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Find gene column
  # ----------------------------------------------------------
  
  gene_candidates <-
    
    grep(
      "gene.*id|geneid|gene_id|^gene$|^id$|symbol",
      colnames(dat),
      ignore.case = TRUE,
      value = TRUE
    )
  
  
  if (length(gene_candidates) > 0) {
    
    gene_col <-
      gene_candidates[1]
    
  } else {
    
    gene_col <-
      colnames(dat)[1]
    
  }
  
  
  # ----------------------------------------------------------
  # Find log2FC
  # ----------------------------------------------------------
  
  lfc_candidates <-
    
    grep(
      "log2FoldChange|log2FC|logFC",
      colnames(dat),
      ignore.case = TRUE,
      value = TRUE
    )
  
  
  if (length(lfc_candidates) == 0) {
    
    stop(
      paste0(
        "\nlog2FoldChange column not found in:\n",
        file,
        "\n\nAvailable columns:\n",
        paste(
          colnames(dat),
          collapse = ", "
        )
      )
    )
    
  }
  
  
  lfc_col <-
    lfc_candidates[1]
  
  
  # ----------------------------------------------------------
  # Find adjusted P value
  # ----------------------------------------------------------
  
  padj_candidates <-
    
    grep(
      "^padj$|adjusted.*p|adj.*p|FDR",
      colnames(dat),
      ignore.case = TRUE,
      value = TRUE
    )
  
  
  if (length(padj_candidates) == 0) {
    
    stop(
      paste0(
        "\npadj/FDR column not found in:\n",
        file,
        "\n\nAvailable columns:\n",
        paste(
          colnames(dat),
          collapse = ", "
        )
      )
    )
    
  }
  
  
  padj_col <-
    padj_candidates[1]
  
  
  # ----------------------------------------------------------
  # Extract
  # ----------------------------------------------------------
  
  result <-
    
    tibble(
      
      Gene =
        as.character(
          dat[[gene_col]]
        ),
      
      log2FC =
        suppressWarnings(
          as.numeric(
            dat[[lfc_col]]
          )
        ),
      
      padj =
        suppressWarnings(
          as.numeric(
            dat[[padj_col]]
          )
        )
      
    )
  
  
  # ----------------------------------------------------------
  # Clean gene ID
  # ----------------------------------------------------------
  
  result$Gene <-
    
    str_extract(
      result$Gene,
      "Solyc[0-9]{2}g[0-9]+(?:\\.[0-9]+)?"
    )
  
  
  result <-
    
    result %>%
    
    filter(
      !is.na(Gene)
    ) %>%
    
    distinct(
      Gene,
      .keep_all = TRUE
    )
  
  
  # ----------------------------------------------------------
  # DEG classification
  # ----------------------------------------------------------
  
  result <-
    
    result %>%
    
    mutate(
      
      DEG =
        
        !is.na(padj) &
        
        padj <= PADJ_CUTOFF &
        
        !is.na(log2FC) &
        
        abs(log2FC) >= LFC_CUTOFF,
      
      
      Direction =
        
        case_when(
          
          DEG &
            log2FC > 0 ~
            "UP",
          
          DEG &
            log2FC < 0 ~
            "DOWN",
          
          TRUE ~
            "NS"
          
        )
      
    )
  
  
  return(result)
  
}


# ============================================================
# 12. READ ALL COMPARISONS
# ============================================================

deseq_data <- list()


for (i in seq_along(deseq_files)) {
  
  
  dat <-
    
    read_deseq(
      deseq_files[i]
    )
  
  
  comparison_letter <-
    
    names(deseq_files)[i]
  
  
  comparison_name <-
    
    sub(
      "_DESeq2_results$",
      "",
      tools::file_path_sans_ext(
        basename(
          deseq_files[i]
        )
      )
    )
  
  
  deseq_data[[i]] <-
    
    list(
      
      letter =
        comparison_letter,
      
      name =
        comparison_name,
      
      file =
        deseq_files[i],
      
      data =
        dat
      
    )
  
  
  cat(
    comparison_letter,
    " | ",
    comparison_name,
    " | DEGs = ",
    sum(
      dat$DEG,
      na.rm = TRUE
    ),
    " | UP = ",
    sum(
      dat$Direction == "UP"
    ),
    " | DOWN = ",
    sum(
      dat$Direction == "DOWN"
    ),
    "\n",
    sep = ""
  )
  
}


# ============================================================
# 13. COMPARISON KEY
# ============================================================

comparison_key <-
  
  tibble(
    
    Letter =
      sapply(
        deseq_data,
        function(x)
          x$letter
      ),
    
    Comparison =
      sapply(
        deseq_data,
        function(x)
          x$name
      )
    
  )


write.csv(
  
  comparison_key,
  
  file.path(
    output_dir,
    "Comparison_Key.csv"
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 14. TOPGO FUNCTION
# ============================================================

run_topGO <- function(
    
  background,
  
  foreground,
  
  gene2go,
  
  ontology,
  
  comparison
  
) {
  
  
  # ----------------------------------------------------------
  # Match genes to GO annotation
  # ----------------------------------------------------------
  
  background <-
    
    intersect(
      background,
      names(gene2go)
    )
  
  
  foreground <-
    
    intersect(
      foreground,
      background
    )
  
  
  if (
    length(background) < 10 ||
    length(foreground) < 2
  ) {
    
    return(NULL)
    
  }
  
  
  # ----------------------------------------------------------
  # Gene vector
  # ----------------------------------------------------------
  
  geneList <-
    
    factor(
      
      as.integer(
        background %in%
          foreground
      )
      
    )
  
  
  names(geneList) <-
    background
  
  
  # ----------------------------------------------------------
  # topGO object
  # ----------------------------------------------------------
  
  GOdata <-
    
    new(
      
      "topGOdata",
      
      ontology =
        ontology,
      
      allGenes =
        geneList,
      
      geneSel =
        function(x)
          x == 1,
      
      annot =
        annFUN.gene2GO,
      
      gene2GO =
        gene2go,
      
      nodeSize =
        5,
      
      description =
        comparison
      
    )
  
  
  # ----------------------------------------------------------
  # Run Fisher test
  # ----------------------------------------------------------
  
  classic_result <-
    
    runTest(
      
      GOdata,
      
      algorithm =
        "classic",
      
      statistic =
        "fisher"
      
    )
  
  
  weight_result <-
    
    runTest(
      
      GOdata,
      
      algorithm =
        "weight01",
      
      statistic =
        "fisher"
      
    )
  
  
  # ----------------------------------------------------------
  # Extract all GO terms
  # ----------------------------------------------------------
  
  res <-
    
    GenTable(
      
      GOdata,
      
      classicFisher =
        classic_result,
      
      weightFisher =
        weight_result,
      
      orderBy =
        "weightFisher",
      
      topNodes =
        length(
          usedGO(GOdata)
        )
      
    )
  
  
  if (
    is.null(res) ||
    nrow(res) == 0
  ) {
    
    return(NULL)
    
  }
  
  
  # ----------------------------------------------------------
  # Numeric p values
  # ----------------------------------------------------------
  
  res$classic_p <-
    
    suppressWarnings(
      as.numeric(
        res$classicFisher
      )
    )
  
  
  res$weight_p <-
    
    suppressWarnings(
      as.numeric(
        res$weightFisher
      )
    )
  
  
  # ----------------------------------------------------------
  # FDR
  # ----------------------------------------------------------
  
  res$FDR <-
    
    p.adjust(
      res$weight_p,
      method = "BH"
    )
  
  
  # ----------------------------------------------------------
  # DEG count
  # ----------------------------------------------------------
  
  res$DEG_Count <-
    
    sapply(
      
      res$GO.ID,
      
      function(go) {
        
        genes <-
          
          genesInTerm(
            GOdata,
            go
          )[[1]]
        
        
        sum(
          genes %in%
            foreground
        )
        
      }
      
    )
  
  
  # ----------------------------------------------------------
  # Background count
  # ----------------------------------------------------------
  
  res$Background_Count <-
    
    sapply(
      
      res$GO.ID,
      
      function(go) {
        
        genes <-
          
          genesInTerm(
            GOdata,
            go
          )[[1]]
        
        
        sum(
          genes %in%
            background
        )
        
      }
      
    )
  
  
  # ----------------------------------------------------------
  # DEG percentage
  # ----------------------------------------------------------
  
  res$DEG_Percent <-
    
    100 *
    
    res$DEG_Count /
    
    length(foreground)
  
  
  # ----------------------------------------------------------
  # -log10 FDR
  # ----------------------------------------------------------
  
  res$negLog10FDR <-
    
    -log10(
      pmax(
        res$FDR,
        1e-300
      )
    )
  
  
  # ----------------------------------------------------------
  # Final table
  # ----------------------------------------------------------
  
  result <-
    
    tibble(
      
      GO_ID =
        res$GO.ID,
      
      GO_name =
        res$Term,
      
      Ontology =
        ontology,
      
      DEG_Count =
        res$DEG_Count,
      
      Background_Count =
        res$Background_Count,
      
      DEG_Percent =
        res$DEG_Percent,
      
      classic_p =
        res$classic_p,
      
      weight_p =
        res$weight_p,
      
      FDR =
        res$FDR,
      
      negLog10FDR =
        res$negLog10FDR
      
    )
  
  
  result
  
}


# ============================================================
# 15. RUN GO ENRICHMENT
# ============================================================

all_results <- list()

enrichment_denominators <- list()   # per Comparison x Ontology totals, used for Fold_Enrichment


for (i in seq_along(deseq_data)) {
  
  
  comp <-
    deseq_data[[i]]
  
  
  dat <-
    comp$data
  
  
  cat(
    "\n========================================\n"
  )
  
  
  cat(
    comp$letter,
    ": ",
    comp$name,
    "\n",
    sep = ""
  )
  
  
  cat(
    "========================================\n"
  )
  
  
  # ----------------------------------------------------------
  # Background
  # ----------------------------------------------------------
  
  background <-
    
    intersect(
      dat$Gene,
      names(gene2GO)
    )
  
  
  # ----------------------------------------------------------
  # UP genes
  # ----------------------------------------------------------
  
  up_genes <-
    
    dat$Gene[
      dat$Direction == "UP"
    ]
  
  
  # ----------------------------------------------------------
  # DOWN genes
  # ----------------------------------------------------------
  
  down_genes <-
    
    dat$Gene[
      dat$Direction == "DOWN"
    ]
  
  
  # ----------------------------------------------------------
  # Save gene lists
  # ----------------------------------------------------------
  
  writeLines(
    
    background,
    
    file.path(
      
      individual_dir,
      
      paste0(
        comp$letter,
        "_background.txt"
      )
      
    )
    
  )
  
  
  writeLines(
    
    up_genes,
    
    file.path(
      
      individual_dir,
      
      paste0(
        comp$letter,
        "_UP.txt"
      )
      
    )
    
  )
  
  
  writeLines(
    
    down_genes,
    
    file.path(
      
      individual_dir,
      
      paste0(
        comp$letter,
        "_DOWN.txt"
      )
      
    )
    
  )
  
  
  # ----------------------------------------------------------
  # Ontology mapping
  # ----------------------------------------------------------
  
  mappings <- list(
    
    BP =
      gene2GO_BP,
    
    CC =
      gene2GO_CC,
    
    MF =
      gene2GO_MF
    
  )
  
  
  # ----------------------------------------------------------
  # Fold-enrichment denominators (per ontology)
  #
  # Calculated exactly the way the KEGG script does it:
  #   Fold_Enrichment = (DEG_Count / Total_DEGs) /
  #                     (Background_Count / Total_Background)
  #
  # Total_Background = number of background genes annotated in
  #   this ontology (BP/CC/MF have different gene universes).
  # Total_DEGs = number of UP+DOWN DEGs (combined, distinct genes)
  #   annotated in this ontology.
  # ----------------------------------------------------------
  
  for (ontology in c("BP", "CC", "MF")) {
    
    ontology_background <-
      intersect(
        background,
        names(mappings[[ontology]])
      )
    
    ontology_deg <-
      intersect(
        union(up_genes, down_genes),
        ontology_background
      )
    
    enrichment_denominators[[
      paste0(comp$letter, "_", ontology)
    ]] <-
      
      tibble(
        Comparison_Letter = comp$letter,
        Ontology          = ontology,
        Total_Background  = length(ontology_background),
        Total_DEGs        = length(ontology_deg)
      )
    
  }
  
  
  # ----------------------------------------------------------
  # UP and DOWN
  # ----------------------------------------------------------
  
  for (
    direction in
    c(
      "UP",
      "DOWN"
    )
  ) {
    
    
    genes <-
      
      if (
        direction == "UP"
      ) {
        
        up_genes
        
      } else {
        
        down_genes
        
      }
    
    
    for (
      ontology in
      c(
        "BP",
        "CC",
        "MF"
      )
    ) {
      
      
      result <-
        
        tryCatch(
          
          run_topGO(
            
            background =
              background,
            
            foreground =
              genes,
            
            gene2go =
              mappings[[ontology]],
            
            ontology =
              ontology,
            
            comparison =
              comp$letter
            
          ),
          
          error =
            function(e) {
              
              message(
                "topGO error in ",
                comp$letter,
                " ",
                direction,
                " ",
                ontology,
                ": ",
                e$message
              )
              
              NULL
              
            }
          
        )
      
      
      if (!is.null(result)) {
        
        
        result$Comparison_Letter <-
          comp$letter
        
        
        result$Comparison_Name <-
          comp$name
        
        
        result$Direction <-
          direction
        
        
        all_results[[
          paste0(
            comp$letter,
            "_",
            direction,
            "_",
            ontology
          )
        ]] <-
          result
        
        
        # ----------------------------------------------------
        # Save individual result
        # ----------------------------------------------------
        
        write.csv(
          
          result,
          
          file.path(
            
            individual_dir,
            
            paste0(
              
              comp$letter,
              "_",
              direction,
              "_",
              ontology,
              "_GO.csv"
              
            )
            
          ),
          
          row.names =
            FALSE
          
        )
        
      }
      
    }
    
  }
  
}


enrichment_denominators_df <-
  bind_rows(enrichment_denominators)


# ============================================================
# 16. COMBINE RESULTS
# ============================================================

if (
  length(all_results) == 0
) {
  
  stop(
    "\nNo GO enrichment results were generated."
  )
  
}


go_results <-
  
  bind_rows(
    all_results
  )


# ============================================================
# 17. SAVE ALL GO RESULTS
# ============================================================

write.csv(
  
  go_results,
  
  file.path(
    
    output_dir,
    
    "GO_All_Results.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 18. SIGNIFICANT GO RESULTS
# ============================================================

significant_go <-
  
  go_results %>%
  
  filter(
    !is.na(FDR),
    FDR <= 0.05
  )


write.csv(
  
  significant_go,
  
  file.path(
    
    output_dir,
    
    "GO_Significant_FDR_0.05.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 19. FUNCTIONAL ANNOTATION
# ============================================================
#
# TOP 10 TERMS PER BP / CC / MF
#
# Ranking:
# total number of DEGs associated with each GO term.
#
# ------------------------------------------------------------
# IMPORTANT:
#
# UP and DOWN are combined into the same vertical bar.
#
# The total bar height = total genes.
#
# The lower part = DOWN
# The upper part = UP
#
# Label = UP/DOWN
#
# ============================================================


annotation_data_list <- list()


for (i in seq_along(deseq_data)) {
  
  
  comp <-
    deseq_data[[i]]
  
  
  dat <-
    comp$data
  
  
  deg <-
    
    dat %>%
    
    filter(
      Direction %in%
        c(
          "UP",
          "DOWN"
        )
    )
  
  
  if (nrow(deg) == 0) {
    next
  }
  
  
  for (j in seq_len(nrow(deg))) {
    
    
    gene <-
      deg$Gene[j]
    
    
    if (
      !gene %in%
      names(gene2GO)
    ) {
      next
    }
    
    
    terms <-
      gene2GO[[gene]]
    
    
    if (length(terms) == 0) {
      next
    }
    
    
    annotation_data_list[[
      length(annotation_data_list) + 1
    ]] <-
      
      tibble(
        
        Comparison =
          comp$letter,
        
        Direction =
          deg$Direction[j],
        
        Gene =
          gene,
        
        GO_ID =
          terms
        
      )
    
  }
  
}


annotation_gene_data <-
  
  bind_rows(
    annotation_data_list
  )


if (
  nrow(annotation_gene_data) == 0
) {
  
  stop(
    "No GO-annotated DEGs were found for functional annotation."
  )
  
}


# ============================================================
# 20. ADD GO INFORMATION
# ============================================================

annotation_gene_data <-
  
  annotation_gene_data %>%
  
  left_join(
    go_info,
    by = "GO_ID"
  ) %>%
  
  filter(
    !is.na(Ontology)
  )


# ============================================================
# 21. COUNT UP/DOWN GENES
# ============================================================

annotation_counts <-
  
  annotation_gene_data %>%
  
  group_by(
    
    Comparison,
    
    Direction,
    
    Ontology,
    
    GO_ID,
    
    GO_name
    
  ) %>%
  
  summarise(
    
    Gene_Count =
      n_distinct(Gene),
    
    .groups =
      "drop"
    
  )


# ============================================================
# 22. TOP 10 TERMS PER ONTOLOGY
# ============================================================

annotation_ranking <-
  
  annotation_counts %>%
  
  group_by(
    
    Ontology,
    
    GO_ID,
    
    GO_name
    
  ) %>%
  
  summarise(
    
    Total_Genes =
      sum(
        Gene_Count,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  group_by(
    Ontology
  ) %>%
  
  arrange(
    desc(Total_Genes)
  ) %>%
  
  slice_head(
    n =
      TOP_ANNOTATION_TERMS
  ) %>%
  
  ungroup()


write.csv(
  
  annotation_ranking,
  
  file.path(
    
    output_dir,
    
    "Functional_Annotation_Top10.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 23. ADD TOP TERMS TO PLOT DATA
# ============================================================

annotation_plot_data <-
  
  annotation_counts %>%
  
  inner_join(
    
    annotation_ranking %>%
      
      select(
        Ontology,
        GO_ID,
        GO_name
      ),
    
    by = c(
      "Ontology",
      "GO_ID",
      "GO_name"
    )
    
  )


# ============================================================
# 24. CREATE UP/DOWN WIDE DATA
# ============================================================

annotation_wide <-
  
  annotation_plot_data %>%
  
  group_by(
    
    Comparison,
    
    Ontology,
    
    GO_ID,
    
    GO_name
    
  ) %>%
  
  summarise(
    
    UP = sum(
      Gene_Count[
        Direction == "UP"
      ],
      na.rm = TRUE
    ),
    
    DOWN = sum(
      Gene_Count[
        Direction == "DOWN"
      ],
      na.rm = TRUE
    ),
    
    .groups =
      "drop"
    
  ) %>%
  
  mutate(
    
    Total =
      UP + DOWN,
    
    Label =
      paste0(
        UP,
        "/",
        DOWN
      )
    
  )


# ============================================================
# 25. ORDER GO TERMS
# ============================================================

annotation_term_order <-
  
  annotation_ranking %>%
  
  group_by(
    Ontology
  ) %>%
  
  arrange(
    Total_Genes
  ) %>%
  
  mutate(
    
    GO_Label =
      make_go_label(
        GO_name,
        GO_ID,
        Ontology,
        colored = FALSE
      )
    
  ) %>%
  
  ungroup()


annotation_wide <-
  
  annotation_wide %>%
  
  left_join(
    
    annotation_term_order %>%
      
      select(
        Ontology,
        GO_ID,
        GO_Label
      ),
    
    by = c(
      "Ontology",
      "GO_ID"
    )
    
  )


# ============================================================
# 26. LONG FORMAT FOR STACKED BARS
# ============================================================

annotation_long <-
  
  annotation_wide %>%
  
  select(
    
    Comparison,
    
    Ontology,
    
    GO_ID,
    
    GO_name,
    
    GO_Label,
    
    UP,
    
    DOWN
    
  ) %>%
  
  pivot_longer(
    
    cols =
      c(
        DOWN,
        UP
      ),
    
    names_to =
      "Direction",
    
    values_to =
      "Gene_Count"
    
  )


annotation_long$Direction <-
  
  factor(
    
    annotation_long$Direction,
    
    levels =
      c(
        "DOWN",
        "UP"
      )
    
  )


annotation_long$Comparison <-
  
  factor(
    
    annotation_long$Comparison,
    
    levels =
      comparison_key$Letter
    
  )


annotation_long$Ontology <-
  
  factor(
    
    annotation_long$Ontology,
    
    levels =
      ONTOLOGY_LEVELS
    
  )


# ============================================================
# 27. FUNCTIONAL ANNOTATION PLOT
# ============================================================
#
# VERTICAL BARS
#
# Rows:
#   comparisons
#
# Columns:
#   BP
#   CC
#   MF
#
# X axis:
#   GO terms
#
# Y axis:
#   number of genes
#
# ============================================================

functional_annotation_plot <-
  
  ggplot(
    
    annotation_long,
    
    aes(
      
      x =
        GO_Label,
      
      y =
        Gene_Count,
      
      fill =
        Direction
      
    )
    
  ) +
  
  geom_col(
    
    width =
      0.75
    
  ) +
  
  facet_grid(
    
    Comparison ~ Ontology,
    
    scales =
      "free_x",
    
    space =
      "free_x",
    
    labeller =
      labeller(
        Ontology =
          ONTOLOGY_LABELS
      )
    
  ) +
  
  # ----------------------------------------------------------
# Number labels
# ----------------------------------------------------------

geom_text(
  
  data =
    
    annotation_wide,
  
  aes(
    
    x =
      GO_Label,
    
    y =
      Total,
    
    label =
      Label
    
  ),
  
  inherit.aes =
    FALSE,
  
  angle =
    90,
  
  hjust =
    -0.2,
  
  vjust =
    0.5,
  
  size =
    FONT_NUMBER,
  
  fontface =
    "bold"
  
) +
  
  scale_fill_manual(
    
    values = c(
      
      DOWN =
        "steelblue",
      
      UP =
        "red"
      
    ),
    
    labels = c(
      
      DOWN =
        "Down",
      
      UP =
        "Up"
      
    )
    
  ) +
  
  scale_y_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0,
            0.18
          )
      )
    
  ) +
  
  labs(
    
    x =
      NULL,
    
    y =
      "Number of genes",
    
    fill =
      NULL
    
  ) +
  
  theme_bw(
    
    base_size =
      FONT_BASE
    
  ) +
  
  theme(
    
    legend.position =
      "top",
    
    legend.text =
      element_text(
        size =
          FONT_LEGEND
      ),
    
    strip.text =
      element_text(
        size =
          FONT_STRIP,
        face =
          "bold"
      ),
    
    strip.background =
      element_rect(
        fill =
          "white"
      ),
    
    panel.grid =
      element_blank(),
    
    panel.background =
      element_rect(
        fill =
          "grey95"
      ),
    
    axis.text.x =
      element_text(
        
        angle =
          90,
        
        hjust =
          1,
        
        vjust =
          0.5,
        
        size =
          FONT_GO,
        
        face =
          "bold",
        
        colour =
          "black"
        
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.title.y =
      element_text(
        size =
          FONT_AXIS_Y,
        
        face =
          "bold"
      ),
    
    panel.spacing.x =
      unit(
        0.35,
        "lines"
      ),
    
    panel.spacing.y =
      unit(
        0.9,
        "lines"
      )
    
  )


# ============================================================
# 28. SAVE FUNCTIONAL ANNOTATION FIGURE
# ============================================================

ggsave(
  
  filename =
    
    file.path(
      
      figure_dir,
      
      "Figure_Functional_Annotation_Top10.tiff"
      
    ),
  
  plot =
    
    functional_annotation_plot,
  
  width =
    ANNOTATION_WIDTH,
  
  height =
    ANNOTATION_HEIGHT,
  
  units =
    "in",
  
  dpi =
    FIG_DPI,
  
  compression =
    "lzw"
  
)


ggsave(
  
  filename =
    
    file.path(
      
      figure_dir,
      
      "Figure_Functional_Annotation_Top10.pdf"
      
    ),
  
  plot =
    
    functional_annotation_plot,
  
  width =
    ANNOTATION_WIDTH,
  
  height =
    ANNOTATION_HEIGHT,
  
  units =
    "in"
  
)


# ============================================================
# 29. FUNCTIONAL ENRICHMENT
# ============================================================
#
# go_results has one row per Comparison x Direction (UP/DOWN)
# x Ontology x GO term. For this figure UP and DOWN are
# combined into a single "number of genes" value per
# Comparison x GO term, and the more significant (smaller)
# FDR of the two directions is kept.
#
# TOP_ENRICHMENT_TERMS_PER_CATEGORY GO terms are then kept
# per ontology (MF / CC / BP), ranked by total gene count
# summed across all comparisons.
#
# ============================================================

go_results_combined <-
  
  go_results %>%
  
  group_by(
    
    Comparison_Letter,
    
    Comparison_Name,
    
    Ontology,
    
    GO_ID,
    
    GO_name
    
  ) %>%
  
  summarise(
    
    DEG_Count =
      sum(
        DEG_Count,
        na.rm = TRUE
      ),
    
    Background_Count =
      max(
        Background_Count,
        na.rm = TRUE
      ),
    
    FDR =
      min(
        FDR,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  left_join(
    
    enrichment_denominators_df,
    
    by = c(
      "Comparison_Letter",
      "Ontology"
    )
    
  ) %>%
  
  mutate(
    
    negLog10FDR =
      -log10(
        pmax(
          FDR,
          1e-300
        )
      ),
    
    # Fold enrichment, calculated the same way as the KEGG script:
    # (DEG_Count / Total_DEGs) / (Background_Count / Total_Background)
    Fold_Enrichment =
      (DEG_Count / Total_DEGs) /
      (Background_Count / Total_Background)
    
  )


top30_enrichment <-
  
  go_results_combined %>%
  
  group_by(
    
    GO_ID,
    
    GO_name,
    
    Ontology
    
  ) %>%
  
  summarise(
    
    Total_DEG =
      sum(
        DEG_Count,
        na.rm = TRUE
      ),
    
    Best_FDR =
      min(
        FDR,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  group_by(
    Ontology
  ) %>%
  
  arrange(
    
    desc(Total_DEG),
    
    Best_FDR,
    
    .by_group =
      TRUE
    
  ) %>%
  
  slice_head(
    
    n =
      TOP_ENRICHMENT_TERMS_PER_CATEGORY
    
  ) %>%
  
  ungroup()


write.csv(
  
  top30_enrichment,
  
  file.path(
    
    output_dir,
    
    "Functional_Enrichment_Top30.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 30. PREPARE ENRICHMENT PLOT DATA
# ============================================================

enrichment_plot_data <-
  
  go_results_combined %>%
  
  inner_join(
    
    top30_enrichment %>%
      
      select(
        GO_ID,
        Ontology
      ),
    
    by = c(
      "GO_ID",
      "Ontology"
    )
    
  ) %>%
  
  mutate(
    
    GO_Label =
      make_go_label(
        GO_name,
        GO_ID,
        Ontology,
        colored = TRUE
      )
    
  )


# ============================================================
# 31. ORDER ENRICHMENT TERMS
# ============================================================
#
# Terms are grouped by category (MF, then CC, then BP, so
# MF ends up at the top of the horizontal bar chart), and
# ranked by total gene count within each category.
# ============================================================

enrichment_order <-
  
  top30_enrichment %>%
  
  mutate(
    
    Ontology =
      factor(
        Ontology,
        levels =
          rev(ONTOLOGY_LEVELS)
      )
    
  ) %>%
  
  arrange(
    
    Ontology,
    
    Total_DEG
    
  ) %>%
  
  mutate(
    
    GO_Label =
      make_go_label(
        GO_name,
        GO_ID,
        Ontology,
        colored = TRUE
      )
    
  ) %>%
  
  pull(
    GO_Label
  )


enrichment_plot_data$GO_Label <-
  
  factor(
    
    enrichment_plot_data$GO_Label,
    
    levels =
      unique(
        enrichment_order
      )
    
  )


enrichment_plot_data$Comparison <-
  
  factor(
    
    enrichment_plot_data$Comparison_Letter,
    
    levels =
      comparison_key$Letter
    
  )


# ============================================================
# 32. FUNCTIONAL ENRICHMENT PLOT
# ============================================================
#
# HORIZONTAL BARS
#
# Y (left side):
#   GO terms, top 10 per category, grouped MF / CC / BP
#
# X:
#   number of genes (UP + DOWN combined)
#
# Columns:
#   8 comparisons (a-h)
#
# Color:
#   Fold enrichment (matches the KEGG functional enrichment figure)
#
# Number:
#   gene count, printed at the end of each bar
#
# Styling below (fonts, spacing, figure size, color scale) is
# matched directly to Figure 2 in 6_2_KEGG_Analysis.R.
# ============================================================

functional_enrichment_plot <-
  
  ggplot(
    
    enrichment_plot_data,
    
    aes(
      
      x =
        DEG_Count,
      
      y =
        GO_Label,
      
      fill =
        Fold_Enrichment
      
    )
    
  ) +
  
  geom_col(
    
    width =
      0.8
    
  ) +
  
  # ----------------------------------------------------------
# Gene count at the end of the bar
# ----------------------------------------------------------

geom_text(
  
  aes(
    
    label =
      DEG_Count
    
  ),
  
  hjust =
    -0.2,
  
  size =
    FONT_NUMBER_ENRICH,
  
  color =
    "black"
  
) +
  
  facet_grid(
    
    . ~ Comparison,
    
    scales =
      "free_x"
    
  ) +
  
  scale_fill_viridis_c(
    
    option =
      "plasma"
    
  ) +
  
  scale_x_continuous(
    
    expand =
      expansion(
        mult =
          c(
            0,
            0.3
          )
      )
    
  ) +
  
  scale_y_discrete(
    
    expand =
      expansion(
        add =
          0.4
      )
    
  ) +
  
  labs(
    
    x =
      "Number of genes",
    
    y =
      NULL,
    
    fill =
      "Fold\nEnrichment"
    
  ) +
  
  theme_bw(
    
    base_size =
      FONT_BASE_ENRICH
    
  ) +
  
  theme(
    
    strip.text =
      element_text(
        
        size =
          FONT_STRIP_ENRICH,
        
        face =
          "bold"
        
      ),
    
    strip.background =
      element_rect(
        fill =
          "white"
      ),
    
    panel.grid =
      element_blank(),
    
    panel.background =
      element_rect(
        fill =
          "grey95"
      ),
    
    axis.text.x =
      element_text(
        size =
          FONT_AXIS_ENRICH - 1
      ),
    
    axis.ticks.x =
      element_line(),
    
    axis.text.y =
      ggtext::element_markdown(
        size =
          FONT_GO_ENRICH,
        
        face =
          "bold"
      ),
    
    axis.title.x =
      element_text(
        size =
          FONT_STRIP_ENRICH,
        
        face =
          "bold"
      ),
    
    legend.text =
      element_text(
        size =
          FONT_LEGEND_ENRICH
      ),
    
    legend.title =
      element_text(
        size =
          FONT_LEGEND_ENRICH,
        
        face =
          "bold"
      ),
    
    panel.spacing.x =
      unit(
        1.2,
        "lines"
      )
    
  )


# ============================================================
# 33. SAVE FUNCTIONAL ENRICHMENT FIGURE
# ============================================================

ggsave(
  
  filename =
    
    file.path(
      
      figure_dir,
      
      "Figure_Functional_Enrichment_Top30.tiff"
      
    ),
  
  plot =
    
    functional_enrichment_plot,
  
  width =
    ENRICHMENT_WIDTH,
  
  height =
    ENRICHMENT_HEIGHT,
  
  units =
    "in",
  
  dpi =
    FIG_DPI,
  
  compression =
    "lzw"
  
)


ggsave(
  
  filename =
    
    file.path(
      
      figure_dir,
      
      "Figure_Functional_Enrichment_Top30.pdf"
      
    ),
  
  plot =
    
    functional_enrichment_plot,
  
  width =
    ENRICHMENT_WIDTH,
  
  height =
    ENRICHMENT_HEIGHT,
  
  units =
    "in"
  
)


# ============================================================
# 34. SAVE FIGURE DATA
# ============================================================

write.csv(
  
  annotation_wide,
  
  file.path(
    
    output_dir,
    
    "Functional_Annotation_Figure_Data.csv"
    
  ),
  
  row.names =
    FALSE
  
)


write.csv(
  
  enrichment_plot_data,
  
  file.path(
    
    output_dir,
    
    "Functional_Enrichment_Figure_Data.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 35. DEG SUMMARY
# ============================================================

deg_summary <-
  
  map_dfr(
    
    deseq_data,
    
    function(x) {
      
      dat <-
        x$data
      
      
      tibble(
        
        Comparison =
          x$name,
        
        Letter =
          x$letter,
        
        Total_Genes =
          nrow(dat),
        
        DEGs =
          sum(
            dat$DEG,
            na.rm = TRUE
          ),
        
        UP =
          sum(
            dat$Direction ==
              "UP"
          ),
        
        DOWN =
          sum(
            dat$Direction ==
              "DOWN"
          )
        
      )
      
    }
    
  )


write.csv(
  
  deg_summary,
  
  file.path(
    
    output_dir,
    
    "DEG_Summary_All_Comparisons.csv"
    
  ),
  
  row.names =
    FALSE
  
)


# ============================================================
# 36. SAVE SESSION INFORMATION
# ============================================================

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  file.path(
    
    output_dir,
    
    "sessionInfo.txt"
    
  )
  
)


# ============================================================
# 37. FINAL MESSAGE
# ============================================================

cat(
  "\n\n====================================================\n"
)

cat(
  "GO ANALYSIS COMPLETED\n"
)

cat(
  "====================================================\n\n"
)

cat(
  "Comparisons:",
  length(deseq_data),
  "\n"
)

cat(
  "DEG cutoff: padj <=",
  PADJ_CUTOFF,
  "\n"
)

cat(
  "LFC cutoff: |log2FC| >=",
  LFC_CUTOFF,
  "\n"
)

cat(
  "Functional annotation:",
  TOP_ANNOTATION_TERMS,
  "terms per BP/CC/MF\n"
)

cat(
  "Functional enrichment:",
  TOP_ENRICHMENT_TERMS,
  "terms\n"
)

cat(
  "Figure resolution:",
  FIG_DPI,
  "dpi\n\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n\n"
)

cat(
  "Figures:\n",
  figure_dir,
  "\n\n"
)

cat(
  "Comparison key:\n"
)

print(
  comparison_key
)

cat(
  "\n====================================================\n"
)