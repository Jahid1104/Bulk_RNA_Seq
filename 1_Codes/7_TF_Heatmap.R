# ==========================================================
# 7_TF_Heatmap.R by Md Jahid Hasan Jone
# Heatmap of selected genes (log2FC) across selected comparisons
# ==========================================================
# What this script does:
#   1. Reads a gene list csv (columns "Gene Name" and "Gene ID", e.g. transcription factors)
#   2. Reads the DESeq2 results csv of each comparison listed in "selected_comparisons" (section 3)
#   3. Pulls log2FoldChange (and padj, for significance stars) for the listed genes
#      out of each comparison's results table
#   4. Builds one matrix: rows = GeneName_GeneID, columns = comparison
#   5. Draws a pheatmap colored by log2FoldChange, with significance stars
#      (* < 0.05, ** < 0.01, *** < 0.001) printed in each cell
#   6. Saves TF_Heatmap.tiff (1000 dpi, width 6.27 in, height 8 in)
#
# Row order = order of the genes in the gene list file
# Column order = the comparisons, in the order listed in section 3
# Neither axis is clustered, so the order stays exactly as specified.
#
# ----------------------------------------------------------
# HOW TO EDIT THE SCRIPT
# ----------------------------------------------------------
#   Paths (section 2): replace /.../.../ with the path to your project folder. The gene list
#     can also be picked in a window with file.choose() (see the Option 2 lines).
#   Gene list: a csv with the columns "Gene Name" and "Gene ID". The Gene ID must match the
#     gene IDs in the DESeq2 results tables (the version suffix is ignored if there is no exact match).
#   Comparisons (section 3): list the comparison names in "selected_comparisons". Each name must
#     match "name" in the comparisons list of 5_DESeq2_Analysis.R exactly, because the script
#     reads <name>_DESeq2_results.csv. Column letters a, b, c ... follow the order of the list.
#       Add a comparison:    add a line with its name
#       Remove a comparison: delete its line or put # in front of it
#       Every line needs a comma at the end except the last one in the list.
#   Figure size (section 7): width and height are set for about 25 genes and 8 comparisons.
#     Increase height for more genes and width for more comparisons.
# ==========================================================




# -----------------------------
# 0. Packages
# -----------------------------
cran_pkgs <- c("pheatmap", "RColorBrewer", "extrafont")
to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(pheatmap)
library(RColorBrewer)
library(extrafont)




# -----------------------------
# 1. Font (Times New Roman, gene names italic)
# -----------------------------
loadfonts(device = if (.Platform$OS.type == "windows") "win" else "pdf", quiet = TRUE)   # Windows. Use device = "pdf" on Mac/Linux
FONT <- "Times New Roman"                 # <- EDIT: use "serif" here if this font isn't found




# -----------------------------
# 2. Paths  <- EDIT these
# -----------------------------
# Folder with the "<comparison>_DESeq2_results.csv" files from step 5
deseq2_dir <- "/.../.../Bulk_RNA_seq/5_Results/5_DESeq2/Tables/Comparisons"

# Gene list csv (columns "Gene Name" and "Gene ID")
# There are two ways to give the gene list. Use ONE of them and comment out the other:
#   Option 1 (default): direct path. Replace /.../.../ with the path to your project folder.
#   Option 2: file.choose() opens a window to pick the file (needs an interactive R session such
#             as RStudio). To use it, remove the # from the "message" and "file.choose()" lines
#             and put a # in front of the direct-path line.

# Option 1: direct path
gene_file   <- "/.../.../Bulk_RNA_seq/5_Results/7_TF/TF_Genes.csv"

# Option 2: choose the file in a window
#message("Choose the GENE LIST csv file (columns: Gene Name, Gene ID)")
#gene_file <- file.choose()

output_dir  <- "/.../.../Bulk_RNA_seq/5_Results/7_TF"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05   # for the significance stars below




# -----------------------------
# 3. The selected comparisons  <- EDIT if this list changes
# -----------------------------
# Must match the comparison "name" values used in the DESeq2 script
# (and therefore the "<name>_DESeq2_results.csv" file names).
selected_comparisons <- c(
  "CLN1466EA_Flower_24h vs NC123S_Flower_24h",
  "CLN1466EA_Flower_72h vs NC123S_Flower_72h",
  "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h",
  "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h",
  "CLN1466EA_Flower vs NC123S_Flower",
  "CLN1466EA_Leaf vs NC123S_Leaf",
  "CLN1466EA vs NC123S",
  "Flower vs Leaf"
)

# Column labels shown on the heatmap (a, b, c ...) -> full comparison name, saved as a legend csv
comparison_labels <- setNames(letters[seq_along(selected_comparisons)], selected_comparisons)




# -----------------------------
# 4. Load gene list
# -----------------------------
gene_list <- read.csv(gene_file, stringsAsFactors = FALSE, check.names = FALSE)
gene_list <- gene_list[, c("Gene Name", "Gene ID")]
names(gene_list) <- c("GeneName", "GeneID")
gene_list <- gene_list[!is.na(gene_list$GeneID), ]   # drop any blank trailing rows

row_label <- paste0(gene_list$GeneName, "_", gene_list$GeneID)




# -----------------------------
# 5. Load the DESeq2 result tables and pull out the listed genes
# -----------------------------
# Some count matrices drop the ITAG version suffix (e.g. "Solyc01g060400" instead
# of "Solyc01g060400.1"), so if an exact Gene ID match fails, fall back to
# matching with the version suffix stripped from both sides.
strip_version <- function(x) sub("\\.\\d+$", "", x)

lfc_mat  <- matrix(NA_real_, nrow = nrow(gene_list), ncol = length(selected_comparisons),
                   dimnames = list(row_label, selected_comparisons))
padj_mat <- lfc_mat

for (comp_name in selected_comparisons) {
  res_file <- file.path(deseq2_dir, paste0(comp_name, "_DESeq2_results.csv"))
  if (!file.exists(res_file)) stop("Missing DESeq2 results file: ", res_file)
  
  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  
  match_idx <- match(gene_list$GeneID, res_df$Gene)
  missing   <- is.na(match_idx)
  if (any(missing)) {
    fallback <- match(strip_version(gene_list$GeneID[missing]), strip_version(res_df$Gene))
    match_idx[missing] <- fallback
  }
  
  lfc_mat[, comp_name]  <- res_df$log2FoldChange[match_idx]
  padj_mat[, comp_name] <- res_df$padj[match_idx]
  
  n_found <- sum(!is.na(match_idx))
  message(comp_name, ": found ", n_found, " / ", nrow(gene_list), " genes")
}




# -----------------------------
# 6. Significance stars
# -----------------------------
star_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat))
star_mat[!is.na(padj_mat) & padj_mat < 0.001]                       <- "***"
star_mat[!is.na(padj_mat) & padj_mat >= 0.001 & padj_mat < 0.01]    <- "**"
star_mat[!is.na(padj_mat) & padj_mat >= 0.01  & padj_mat < padj_cutoff] <- "*"




# -----------------------------
# 7. Heatmap
# -----------------------------
# Panel letter -> full comparison name, for reference
write.csv(data.frame(Panel = unname(comparison_labels), Comparison = names(comparison_labels)),
          file.path(output_dir, "TF_Heatmap_Column_Legend.csv"), row.names = FALSE)

# Swap in the short letter labels for the plot itself
colnames(lfc_mat)  <- comparison_labels[colnames(lfc_mat)]
colnames(star_mat) <- comparison_labels[colnames(star_mat)]

italic_labels <- as.expression(lapply(rownames(lfc_mat), function(g) bquote(italic(.(g)))))

# symmetric diverging scale so 0 (no change) sits at white
max_abs <- max(abs(lfc_mat), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)
palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

tiff(file.path(output_dir, "TF_Heatmap.tiff"),   # saved directly into 5_Results/7_TF
     width = 6.27, height = 8, units = "in", res = 1000, compression = "lzw")

pheatmap(
  lfc_mat,
  color            = palette,
  breaks           = breaks,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  display_numbers  = star_mat,
  number_color     = "black",
  na_col           = "grey90",
  labels_row       = italic_labels,
  angle_col        = 0,           # single letters (a, b, c ...), horizontal reads fine and saves space
  fontsize         = 10,
  fontsize_row     = 10,
  fontsize_col     = 10,
  fontsize_number  = 10,
  fontfamily       = FONT,
  border_color     = "grey60",
  legend           = TRUE,
  main             = " "
)

dev.off()

message("\nDone. Saved: ", normalizePath(file.path(output_dir, "TF_Heatmap.tiff")))