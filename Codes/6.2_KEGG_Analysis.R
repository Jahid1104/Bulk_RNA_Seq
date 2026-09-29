# ============================================================================
# KEGG PATHWAY ANALYSIS FOR TOMATO RNA-SEQ (8 COMPARISONS a-h)
#
# WHAT THIS SCRIPT DOES, STEP BY STEP
#   1. Reads the KAAS query.ko.txt file (Gene ID -> K number).
#   2. Downloads the official KEGG pathway hierarchy once, so every pathway
#      can be labelled with a category (Metabolism, Genetic Information
#      Processing, Environmental Information Processing, Cellular Processes,
#      Organismal Systems).
#   3. Reads the 8 DESeq2 result tables (a-h) and attaches K numbers to genes.
#   4. Flags significant DEGs (padj < 0.05, |log2FC| > 1) and splits them
#      into up- and down-regulated.
#   5. Counts how many up/down DEGs map to each KEGG pathway
#      (= "functional annotation").
#   6. Calculates, for every pathway and every comparison, the fold
#      enrichment: (DEG_Count / total DEGs) / (Background_Count / total
#      background genes)  (= "functional enrichment"). A formal statistical
#      enrichKEGG() test is also run and saved to CSV for the record, but the
#      figure itself is built from the direct calculation above so that every
#      one of the requested pathways is always shown, not only the ones that
#      happened to pass a significance cutoff.
#   7. Draws the two summary figures, styled the same way as the Gene
#      Ontology figures (grey panel background, black plain-text labels on
#      the annotation figure, category-coloured labels on the enrichment
#      figure, "Up/Down" or DEG-count numbers printed on every bar).
#
# OUTPUT
#   7_KEGG/
#     01_KO_annotation/        Gene ID <-> K number table
#     02_Annotated_DESeq2/     DESeq2 results + K number + category
#     03_DEG_tables/           Significant / up / down DEG tables per comparison
#     04_KEGG_enrichment/      enrichKEGG() statistical result tables (for the record)
#     Figure1_Functional_Annotation.tiff
#     Figure2_Functional_Enrichment.tiff
#     Comparison_Key.csv       what letter a-h stands for
#
# NOTE ON GENE ID MATCHING
#   ITAG4.0 transcript IDs look like "Solyc10g079470.3.1" (gene.version.mRNA).
#   DESeq2 was run on gene-level counts, so its IDs usually drop the last
#   ".1" (mRNA number) but keep the gene version, e.g. "Solyc10g079470.3".
#   Rather than assume this, the script tests a couple of ways of trimming
#   the KAAS gene ID and automatically keeps whichever matches your DESeq2
#   gene IDs best (see section 6). Check the printed match rate the first
#   time you run this.
#
# WHY FIGURE 1 LABELS ARE PLAIN BLACK BUT FIGURE 2 LABELS ARE COLOURED
#   This mirrors the Gene Ontology script exactly. Rotated, category-coloured
#   markdown text (as would be needed on Figure 1's x-axis) is not reliably
#   rendered by ggtext, so Figure 1 uses plain black rotated text instead
#   (category is still visible from the column it sits in). Figure 2's
#   labels are horizontal, where coloured markdown renders reliably.
# ============================================================================


############################################################
# 1. LOAD PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

needed_cran <- c("dplyr", "tidyr", "readr", "stringr", "ggplot2", "purrr", "tibble", "ggtext", "jsonlite", "ragg")
for (pkg in needed_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

needed_bioc <- c("clusterProfiler", "KEGGREST")
for (pkg in needed_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(purrr)
library(tibble)
library(ggtext)
library(jsonlite)
library(clusterProfiler)
library(KEGGREST)

# clusterProfiler masks dplyr::select() / dplyr::filter() -- make sure the
# dplyr versions win, since that is what the rest of this script expects.
select <- dplyr::select
filter <- dplyr::filter


############################################################
# 2. USER SETTINGS
############################################################

# ---- input files -----------------------------------------------------------
ko_file   <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/Reference/query.ko.txt"
deseq_dir <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/New_Name/Results/5_DESeq2/Tables/Comparisons"

# ---- output location --------------------------------------------------------
output_dir <- "R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/New_Name/Results/7_KEGG"

# ---- the 8 comparisons, labelled a-h ----------------------------------------
comparison_files <- c(
  a = "CLN1466EA vs NC123S_DESeq2_results.csv",
  b = "Flower vs Leaf_DESeq2_results.csv",
  c = "CLN1466EA_Flower vs NC123S_Flower_DESeq2_results.csv",
  d = "CLN1466EA_Flower_24h vs NC123S_Flower_24h_DESeq2_results.csv",
  e = "CLN1466EA_Flower_72h vs NC123S_Flower_72h_DESeq2_results.csv",
  f = "CLN1466EA_Leaf vs NC123S_Leaf_DESeq2_results.csv",
  g = "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h_DESeq2_results.csv",
  h = "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h_DESeq2_results.csv"
)

# ---- DEG thresholds ----------------------------------------------------------
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1   # |log2FoldChange| > 1  ==  2-fold change

# ---- KEGG enrichment settings (used for the statistical record tables only) --
MIN_GS_SIZE <- 3
MAX_GS_SIZE <- 500

# ---- how many top pathways to show on each figure -----------------------------
TOP_N_ANNOTATION  <- 30   # Figure 1
TOP_N_ENRICHMENT  <- 20   # Figure 2

# ---- the 5 KEGG categories this analysis is restricted to ---------------------
KEEP_CATEGORIES <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Organismal Systems"
)

# ---- colours -------------------------------------------------------------------
UP_COLOR   <- "red"
DOWN_COLOR <- "steelblue"

CATEGORY_COLORS <- c(
  "Metabolism"                            = "#1B9E77",
  "Genetic Information Processing"        = "#D95F02",
  "Environmental Information Processing"  = "#7570B3",
  "Cellular Processes"                    = "#E7298A",
  "Organismal Systems"                    = "#66A61E"
)

# ---- font sizes (content >= 9, categories/titles 12-14) -----------------------
FONT_BASE    <- 9    # overall base size
FONT_AXIS    <- 9    # axis text (numeric ticks, gene counts)
FONT_PATHWAY <- 10   # pathway name labels (bold, larger than the default axis text)
FONT_STRIP   <- 12   # comparison a-h strip labels
FONT_LEGEND  <- 9    # legend text
FONT_NUMBER  <- 2.6  # geom_text data-label size (ggplot geom_text uses mm, roughly x2.8 = pt)

# ---- figure sizes (inches) -----------------------------------------------------
FIG1_WIDTH  <- 12
FIG1_HEIGHT <- 18
FIG2_WIDTH  <- 18
FIG2_HEIGHT <- 12
FIG_DPI     <- 1000


############################################################
# 3. CREATE OUTPUT FOLDERS
############################################################

subfolders <- c("01_KO_annotation", "02_Annotated_DESeq2", "03_DEG_tables", "04_KEGG_enrichment")
for (sub in subfolders) {
  dir.create(file.path(output_dir, sub), showWarnings = FALSE, recursive = TRUE)
}


############################################################
# 4. HELPER FUNCTIONS
############################################################

# ---- read the KAAS query.ko.txt file ----
# Some genes have no K number, so lines can have 1 or 2 columns.
# Reading line-by-line and splitting manually handles that safely.
read_kaas_ko <- function(path) {
  
  if (!file.exists(path)) {
    stop("KO number file not found:\n", path)
  }
  
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\t")
  
  tibble(
    Gene_ID  = str_trim(vapply(parts, `[`, character(1), 1)),
    K_Number = str_trim(vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)))
  ) %>%
    mutate(K_Number = str_extract(K_Number, "K[0-9]{5}")) %>%
    filter(Gene_ID != "", !is.na(Gene_ID)) %>%
    distinct(Gene_ID, K_Number)
}

# ---- remove the last `n` dot-separated segments from an ID ----
# e.g. strip_id_segments("Solyc10g079470.3.1", 1) -> "Solyc10g079470.3"
#      strip_id_segments("Solyc10g079470.3.1", 2) -> "Solyc10g079470"
strip_id_segments <- function(x, n) {
  if (n == 0) return(x)
  pattern <- paste0("(\\.[^.]+){", n, "}$")
  sub(pattern, "", x)
}

# ---- figure out which ID-trimming strategy matches your DESeq2 gene IDs ----
# Tries no trimming, trimming 1 segment, and trimming 2 segments, and keeps
# whichever gives the highest overlap with a real DESeq2 gene ID column.
choose_id_match_strategy <- function(kaas_gene_ids, deseq_gene_ids) {
  
  match_rates <- sapply(0:2, function(n) {
    trimmed <- unique(strip_id_segments(kaas_gene_ids, n))
    mean(deseq_gene_ids %in% trimmed)
  })
  
  best_n <- (0:2)[which.max(match_rates)]
  
  cat("\nGene ID matching check (against comparison 'a'):\n")
  for (n in 0:2) {
    cat("  trim", n, "segment(s): ", round(match_rates[n + 1] * 100, 1), "% matched\n", sep = "")
  }
  cat("  -> using trim =", best_n, "\n\n")
  
  if (max(match_rates) < 0.5) {
    warning("Less than 50% of DESeq2 gene IDs matched the KO annotation. ",
            "Double check the ID formats before trusting the results.")
  }
  
  best_n
}

# ---- read a DESeq2 result file and standardise the gene ID column ----
read_deseq_file <- function(path) {
  
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  
  gene_col_candidates <- c("GeneID", "Gene_ID", "gene_id", "Gene", "gene", "X", "ID")
  gene_col <- intersect(gene_col_candidates, colnames(df))
  
  if (length(gene_col) > 0) {
    df <- df %>% rename(Gene_ID = all_of(gene_col[1]))
  } else {
    # no obvious gene ID column -- assume the first column is the gene ID
    colnames(df)[1] <- "Gene_ID"
  }
  
  required <- c("log2FoldChange", "padj")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("File is missing required column(s): ", paste(missing, collapse = ", "),
         "\nFile: ", path)
  }
  
  df
}

# small helper (base R has no built-in %||%)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- download (once) the official KEGG pathway -> category map ----
# Source: KEGG BRITE hierarchy "ko00001" (Metabolism, Genetic Information
# Processing, Environmental Information Processing, Cellular Processes,
# Organismal Systems, Human Diseases, Drug Development).
get_kegg_category_map <- function(cache_file) {
  
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  cat("Downloading KEGG pathway hierarchy (one-time download)...\n")
  brite <- fromJSON("https://rest.kegg.jp/get/br:ko00001/json", simplifyVector = FALSE)
  
  rows <- list()
  
  for (category in brite$children) {
    category_name <- str_remove(category$name, "^[0-9]+\\s+")
    
    for (subcategory in category$children) {
      for (pathway in subcategory$children %||% list()) {
        
        pathway_number <- str_extract(pathway$name, "^[0-9]{5}")
        if (is.na(pathway_number)) next
        
        pathway_name <- pathway$name %>%
          str_remove("^[0-9]{5}\\s+") %>%
          str_remove("\\s*\\[PATH:.*\\]$")
        
        rows[[length(rows) + 1]] <- tibble(
          Pathway_Number = pathway_number,
          Pathway_Name   = pathway_name,
          Category       = category_name
        )
      }
    }
  }
  
  category_map <- bind_rows(rows) %>% distinct(Pathway_Number, .keep_all = TRUE)
  saveRDS(category_map, cache_file)
  category_map
}

# ---- download (once) the K number -> pathway map ----
get_ko_pathway_map <- function(cache_file) {
  
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  
  cat("Downloading K number -> pathway map (one-time download)...\n")
  link <- keggLink("pathway", "ko")
  
  result <- tibble(
    K_Number       = str_remove(names(link), "^ko:"),
    Pathway_Number = str_extract(unname(link), "[0-9]{5}")
  ) %>%
    filter(!is.na(Pathway_Number)) %>%
    distinct()
  
  saveRDS(result, cache_file)
  result
}

# ---- shorten a pathway name to `n_words` words, adding "..." if cut ----
truncate_name <- function(x, n_words = 3) {
  vapply(str_split(x, "\\s+"), function(words) {
    if (length(words) <= n_words) {
      paste(words, collapse = " ")
    } else {
      paste0(paste(words[1:n_words], collapse = " "), " ...")
    }
  }, character(1))
}

# ---- count DISTINCT genes per pathway for a set of DEGs ----
# A DEG that maps to several pathways is counted once per pathway; a DEG
# whose K number is shared by another DEG still counts as its own gene.
count_genes_per_pathway <- function(gene_ids, ko_annotation, ko_pathway_map, category_map) {
  
  tibble(Gene_ID = gene_ids) %>%
    left_join(ko_annotation, by = c("Gene_ID" = "Gene_ID_matched")) %>%
    filter(!is.na(K_Number)) %>%
    left_join(ko_pathway_map, by = "K_Number", relationship = "many-to-many") %>%
    filter(!is.na(Pathway_Number)) %>%
    inner_join(category_map, by = "Pathway_Number") %>%   # keeps only the 5 target categories
    distinct(Gene_ID, Pathway_Number, Pathway_Name, Category) %>%
    count(Pathway_Number, Pathway_Name, Category, name = "Gene_Count")
}


############################################################
# 5. READ KO ANNOTATION AND BUILD KEGG REFERENCE TABLES
############################################################

cat("==================================================\n")
cat("STEP 1: KO annotation and KEGG reference tables\n")
cat("==================================================\n")

ko_annotation <- read_kaas_ko(ko_file)

write.csv(
  ko_annotation,
  file.path(output_dir, "01_KO_annotation", "Gene_KNumber_table.csv"),
  row.names = FALSE
)

cat("Genes with a K number:", sum(!is.na(ko_annotation$K_Number)), "of", nrow(ko_annotation), "\n")

# all K numbers assigned in the reference proteome = enrichment background
kegg_background <- ko_annotation %>% filter(!is.na(K_Number)) %>% pull(K_Number) %>% unique()
cat("Background K numbers:", length(kegg_background), "\n")

# pathway -> category, and K number -> pathway (cached locally so repeat
# runs don't re-download from KEGG every time)
category_map   <- get_kegg_category_map(file.path(output_dir, "kegg_category_map.rds"))
ko_pathway_map <- get_ko_pathway_map(file.path(output_dir, "ko_pathway_map.rds"))

# keep only the 5 categories requested for this analysis
category_map <- category_map %>%
  filter(Category %in% KEEP_CATEGORIES) %>%
  mutate(Category = factor(Category, levels = KEEP_CATEGORIES))


############################################################
# 6. CHOOSE THE GENE ID MATCHING STRATEGY
############################################################

first_deseq <- read_deseq_file(file.path(deseq_dir, comparison_files[["a"]]))
id_trim_n <- choose_id_match_strategy(ko_annotation$Gene_ID, first_deseq$Gene_ID)

ko_annotation <- ko_annotation %>%
  mutate(Gene_ID_matched = strip_id_segments(Gene_ID, id_trim_n)) %>%
  distinct(Gene_ID_matched, K_Number)


############################################################
# 7. PROCESS EACH COMPARISON: ANNOTATE + FLAG DEGs
############################################################

cat("\n==================================================\n")
cat("STEP 2: annotating DESeq2 results and flagging DEGs\n")
cat("==================================================\n")

deg_results <- list()   # significant DEGs (all/up/down) per comparison

for (label in names(comparison_files)) {
  
  cat("\nComparison", label, "-", comparison_files[[label]], "\n")
  
  deseq_path <- file.path(deseq_dir, comparison_files[[label]])
  deseq_df   <- read_deseq_file(deseq_path)
  
  annotated <- deseq_df %>%
    left_join(ko_annotation, by = c("Gene_ID" = "Gene_ID_matched")) %>%
    left_join(ko_pathway_map, by = "K_Number", relationship = "many-to-many") %>%
    left_join(category_map, by = "Pathway_Number")
  
  write.csv(
    annotated,
    file.path(output_dir, "02_Annotated_DESeq2", paste0(label, "_annotated.csv")),
    row.names = FALSE
  )
  
  sig <- deseq_df %>%
    filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_CUTOFF)
  
  up   <- sig %>% filter(log2FoldChange > 0)
  down <- sig %>% filter(log2FoldChange < 0)
  
  write.csv(sig,  file.path(output_dir, "03_DEG_tables", paste0(label, "_all_DEGs.csv")),  row.names = FALSE)
  write.csv(up,   file.path(output_dir, "03_DEG_tables", paste0(label, "_UP_DEGs.csv")),   row.names = FALSE)
  write.csv(down, file.path(output_dir, "03_DEG_tables", paste0(label, "_DOWN_DEGs.csv")), row.names = FALSE)
  
  deg_results[[label]] <- list(all = sig, up = up, down = down)
  
  cat("  significant DEGs:", nrow(sig), " (up:", nrow(up), ", down:", nrow(down), ")\n")
}


############################################################
# 8. DEG COUNTS PER PATHWAY (UP vs DOWN), ALL COMPARISONS
############################################################
# This is the shared building block for both figures: how many up- and
# down-regulated genes from each comparison map to each pathway.

annotation_counts <- map_dfr(names(deg_results), function(label) {
  
  up_counts <- count_genes_per_pathway(deg_results[[label]]$up$Gene_ID, ko_annotation, ko_pathway_map, category_map) %>%
    mutate(Direction = "UP")
  
  down_counts <- count_genes_per_pathway(deg_results[[label]]$down$Gene_ID, ko_annotation, ko_pathway_map, category_map) %>%
    mutate(Direction = "DOWN")
  
  bind_rows(up_counts, down_counts) %>% mutate(Comparison = label)
})


############################################################
# 9. RANK ALL PATHWAYS BY TOTAL DEGs ANNOTATED (ALL COMPARISONS)
############################################################

pathway_ranking <- annotation_counts %>%
  group_by(Pathway_Number, Pathway_Name, Category) %>%
  summarise(Total_Genes = sum(Gene_Count), .groups = "drop") %>%
  arrange(desc(Total_Genes))

write.csv(pathway_ranking, file.path(output_dir, "Pathway_Ranking_AllComparisons.csv"), row.names = FALSE)

# ---- top pathways for each figure (Figure 1 gets more than Figure 2) ----
build_pathway_labels <- function(ranking, n) {
  
  ranking %>%
    slice_head(n = n) %>%
    arrange(Category, desc(Total_Genes)) %>%             # group by category for the axis order
    mutate(
      Short_Name    = truncate_name(Pathway_Name, 3),
      Label_Colored = paste0("<span style='color:", CATEGORY_COLORS[as.character(Category)], "'>", Short_Name, "</span>"),
      Pathway_Number = factor(Pathway_Number, levels = Pathway_Number)   # locks the plotting order
    )
}

top_annotation <- build_pathway_labels(pathway_ranking, TOP_N_ANNOTATION)
top_enrichment <- build_pathway_labels(pathway_ranking, TOP_N_ENRICHMENT)

write.csv(top_annotation, file.path(output_dir, "Top_Pathways_Figure1.csv"), row.names = FALSE)
write.csv(top_enrichment, file.path(output_dir, "Top_Pathways_Figure2.csv"), row.names = FALSE)


############################################################
# 10. FOLD ENRICHMENT (used directly for Figure 2)
############################################################
# Fold enrichment = (DEG_Count / total DEGs) / (Background_Count / total
# background genes), calculated directly so every one of the top pathways
# gets a value for every comparison -- it does not depend on a pathway
# having passed a significance cutoff.

cat("\n==================================================\n")
cat("STEP 3: fold enrichment per comparison\n")
cat("==================================================\n")

# background: how many reference genes map to each pathway, and how many
# reference genes map to our 5 kept categories in total
background_counts <- count_genes_per_pathway(unique(ko_annotation$Gene_ID_matched), ko_annotation, ko_pathway_map, category_map) %>%
  rename(Background_Count = Gene_Count)

total_background_genes <- ko_annotation %>%
  filter(!is.na(K_Number)) %>%
  left_join(ko_pathway_map, by = "K_Number", relationship = "many-to-many") %>%
  filter(!is.na(Pathway_Number)) %>%
  inner_join(category_map, by = "Pathway_Number") %>%
  distinct(Gene_ID_matched) %>%
  nrow()

fold_enrichment_all <- map_dfr(names(deg_results), function(label) {
  
  deg_gene_ids <- deg_results[[label]]$all$Gene_ID
  
  deg_counts <- count_genes_per_pathway(deg_gene_ids, ko_annotation, ko_pathway_map, category_map) %>%
    rename(DEG_Count = Gene_Count)
  
  total_deg_genes <- tibble(Gene_ID = deg_gene_ids) %>%
    left_join(ko_annotation, by = c("Gene_ID" = "Gene_ID_matched")) %>%
    filter(!is.na(K_Number)) %>%
    left_join(ko_pathway_map, by = "K_Number", relationship = "many-to-many") %>%
    filter(!is.na(Pathway_Number)) %>%
    inner_join(category_map, by = "Pathway_Number") %>%
    distinct(Gene_ID) %>%
    nrow()
  
  deg_counts %>%
    left_join(background_counts, by = c("Pathway_Number", "Pathway_Name", "Category")) %>%
    mutate(
      Comparison        = label,
      Total_DEGs        = total_deg_genes,
      Total_Background  = total_background_genes,
      Fold_Enrichment   = (DEG_Count / Total_DEGs) / (Background_Count / Total_Background)
    )
})

write.csv(
  fold_enrichment_all,
  file.path(output_dir, "04_KEGG_enrichment", "Fold_Enrichment_AllComparisons.csv"),
  row.names = FALSE
)


############################################################
# 11. STATISTICAL KEGG ENRICHMENT (enrichKEGG, saved for the record)
############################################################
# Kept as a separate statistical record (p-values, FDR) alongside the
# direct fold-enrichment calculation used for the figure above.

for (label in names(deg_results)) {
  
  deg_k_numbers <- deg_results[[label]]$all %>%
    left_join(ko_annotation, by = c("Gene_ID" = "Gene_ID_matched")) %>%
    filter(!is.na(K_Number)) %>%
    pull(K_Number) %>%
    unique()
  
  if (length(deg_k_numbers) < 2) {
    cat("Comparison", label, ": too few annotated DEGs for enrichKEGG(), skipping.\n")
    next
  }
  
  result <- tryCatch(
    enrichKEGG(
      gene          = deg_k_numbers,
      organism      = "ko",
      keyType       = "kegg",
      universe      = kegg_background,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.20,
      minGSSize     = MIN_GS_SIZE,
      maxGSSize     = MAX_GS_SIZE
    ),
    error = function(e) {
      warning("enrichKEGG failed for comparison ", label, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    cat("Comparison", label, ": no significantly enriched pathways (enrichKEGG).\n")
    next
  }
  
  write.csv(
    as.data.frame(result),
    file.path(output_dir, "04_KEGG_enrichment", paste0(label, "_enrichKEGG_statistical.csv")),
    row.names = FALSE
  )
  
  cat("Comparison", label, ":", nrow(as.data.frame(result)), "significantly enriched pathways.\n")
}


############################################################
# 12. SHARED PLOT STYLING (matches the GO script)
############################################################

comparison_levels <- names(comparison_files)

base_theme <- theme_bw(base_size = FONT_BASE) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "grey95"),   # same shade as the GO figures
    strip.background = element_rect(fill = "white"),
    strip.text       = element_text(size = FONT_STRIP, face = "bold"),
    axis.text        = element_text(size = FONT_AXIS, color = "black"),
    axis.title       = element_text(size = FONT_STRIP, face = "bold"),
    legend.title     = element_text(size = FONT_STRIP - 1, face = "bold"),
    legend.text      = element_text(size = FONT_LEGEND)
  )


############################################################
# 13. FIGURE 1 - FUNCTIONAL ANNOTATION (up/down gene counts)
############################################################

cat("\nBuilding Figure 1: Functional Annotation...\n")

# one row per Comparison x Pathway x Direction (UP/DOWN), summed across genes.
# Pathway/Comparison combos with no genes at all simply have no bar -- the
# x-axis still shows every one of the top pathways because Label_Colored is
# a factor with all TOP_N_ANNOTATION levels fixed in advance.
fig1_data <- annotation_counts %>%
  filter(Pathway_Number %in% levels(top_annotation$Pathway_Number)) %>%
  group_by(Comparison, Pathway_Number, Direction) %>%
  summarise(Gene_Count = sum(Gene_Count), .groups = "drop") %>%
  left_join(top_annotation %>% select(Pathway_Number, Label_Colored), by = "Pathway_Number") %>%
  mutate(
    Comparison    = factor(Comparison, levels = comparison_levels),
    Direction     = factor(Direction, levels = c("DOWN", "UP")),
    Label_Colored = factor(Label_Colored, levels = unique(top_annotation$Label_Colored)),
    Signed_Count  = if_else(Direction == "UP", Gene_Count, -Gene_Count)   # UP = +, DOWN = -
  )

figure1 <- ggplot(fig1_data, aes(x = Label_Colored, y = Signed_Count, fill = Direction)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  geom_col(position = "identity", width = 0.8) +
  geom_text(
    data = fig1_data %>% filter(Gene_Count > 0),
    aes(label = Gene_Count, hjust = if_else(Direction == "UP", -0.3, 1.3)),
    angle = 90, vjust = 0.5, size = FONT_NUMBER, fontface = "bold"
  ) +
  facet_grid(rows = vars(Comparison)) +
  scale_fill_manual(values = c(UP = UP_COLOR, DOWN = DOWN_COLOR), labels = c(UP = "Up", DOWN = "Down")) +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.25))) +
  labs(x = NULL, y = "Number of genes (Up = +, Down = -)", fill = NULL) +
  base_theme +
  theme(
    legend.position = "top",
    axis.text.x     = ggtext::element_markdown(angle = 90, hjust = 1, vjust = 0.5, size = FONT_PATHWAY, face = "bold"),
    panel.spacing   = unit(0.4, "lines")
  )

ggsave(
  file.path(output_dir, "Figure1_Functional_Annotation.tiff"),
  figure1, width = FIG1_WIDTH, height = FIG1_HEIGHT, dpi = FIG_DPI,
  compression = "lzw", device = ragg::agg_tiff
)


############################################################
# 14. FIGURE 2 - FUNCTIONAL ENRICHMENT (fold enrichment)
############################################################

cat("Building Figure 2: Functional Enrichment...\n")

fig2_data <- fold_enrichment_all %>%
  filter(Pathway_Number %in% levels(top_enrichment$Pathway_Number)) %>%
  left_join(top_enrichment %>% select(Pathway_Number, Label_Colored), by = "Pathway_Number") %>%
  mutate(
    Comparison     = factor(Comparison, levels = comparison_levels),
    Pathway_Number = factor(Pathway_Number, levels = rev(levels(top_enrichment$Pathway_Number))),
    Label_Colored  = factor(Label_Colored, levels = rev(unique(top_enrichment$Label_Colored)))
  )

figure2 <- ggplot(fig2_data, aes(x = DEG_Count, y = Label_Colored, fill = Fold_Enrichment)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = DEG_Count), hjust = -0.2, size = FONT_NUMBER, color = "black") +
  facet_grid(cols = vars(Comparison)) +
  scale_fill_viridis_c(option = "plasma") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  scale_y_discrete(expand = expansion(add = 0.4)) +   # tight gap between pathway rows
  labs(x = "Number of DEGs", y = NULL, fill = "Fold\nEnrichment") +
  base_theme +
  theme(
    axis.text.y     = ggtext::element_markdown(size = FONT_PATHWAY, face = "bold"),
    axis.text.x     = element_text(size = FONT_AXIS - 1),
    panel.spacing.x = unit(1.2, "lines")   # wider gap so the DEG-count axis doesn't overlap between panels
  )

ggsave(
  file.path(output_dir, "Figure2_Functional_Enrichment.tiff"),
  figure2, width = FIG2_WIDTH, height = FIG2_HEIGHT, dpi = FIG_DPI,
  compression = "lzw", device = ragg::agg_tiff
)


############################################################
# 15. SAVE COMPARISON KEY AND FINISH
############################################################

write.csv(
  tibble(Letter = names(comparison_files), Comparison = unname(comparison_files)),
  file.path(output_dir, "Comparison_Key.csv"),
  row.names = FALSE
)

cat("\n==================================================\n")
cat("DONE. Output written to:\n", output_dir, "\n")
cat("==================================================\n")