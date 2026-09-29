# ==========================================================
# HSP70_Heatmap.R by Md Jahid Hasan Jone
# Heatmap of 25 selected HSP70 genes (log2FC) across 8 comparisons
# ==========================================================
# What this script does:
#   1. Reads the HSP70 gene list (Gene Name, Gene ID) from the xlsx file
#   2. Reads the DESeq2 results csv for each of the 8 selected comparisons
#   3. Pulls log2FoldChange (and padj, for significance stars) for the 25 genes
#      out of each comparison's results table
#   4. Builds one matrix: rows = GeneName_GeneID, columns = comparison
#   5. Draws a pheatmap colored by log2FoldChange, with significance stars
#      (* < 0.05, ** < 0.01, *** < 0.001) printed in each cell
#   6. Saves HSP70_Heatmap.tiff (1000 dpi, width 6.27 in, height 8 in)
#
# Row order = HSP70-1 ... HSP70-25 (as listed in the gene list file)
# Column order = the 8 comparisons, in the order listed below
# Neither axis is clustered, so the order stays exactly as specified.
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
loadfonts(device = "win", quiet = TRUE)   # Windows. Use device = "pdf" on Mac/Linux
FONT <- "Times New Roman"                 # <- EDIT: use "serif" here if this font isn't found




# -----------------------------
# 2. Paths  <- EDIT these
# -----------------------------
setwd("R:/Md_Jahid_Hasan_Jone/Experiments_and_Data/5_RNA_seq/New_Name")

deseq2_dir  <- "Results/5_DESeq2/Tables/Comparisons"   # folder with the "<comparison>_DESeq2_results.csv" files
hsp_file    <- "Results/7_HSP/HSP_Genes.csv"           # gene list: columns "Gene Name", "Gene ID"

output_dir  <- "Results/7_HSP"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05   # for the significance stars below




# -----------------------------
# 3. The 8 selected comparisons  <- EDIT if this list changes
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

# Column labels shown on the heatmap (a-h) -> full comparison name, saved as a legend csv
comparison_labels <- setNames(letters[seq_along(selected_comparisons)], selected_comparisons)




# -----------------------------
# 4. Load HSP70 gene list
# -----------------------------
hsp_genes <- read.csv(hsp_file, stringsAsFactors = FALSE, check.names = FALSE)
hsp_genes <- hsp_genes[, c("Gene Name", "Gene ID")]
names(hsp_genes) <- c("GeneName", "GeneID")
hsp_genes <- hsp_genes[!is.na(hsp_genes$GeneID), ]   # drop any blank trailing rows

row_label <- paste0(hsp_genes$GeneName, "_", hsp_genes$GeneID)




# -----------------------------
# 5. Load the 8 DESeq2 result tables and pull out the HSP70 genes
# -----------------------------
# Some count matrices drop the ITAG version suffix (e.g. "Solyc01g060400" instead
# of "Solyc01g060400.1"), so if an exact Gene ID match fails, fall back to
# matching with the version suffix stripped from both sides.
strip_version <- function(x) sub("\\.\\d+$", "", x)

lfc_mat  <- matrix(NA_real_, nrow = nrow(hsp_genes), ncol = length(selected_comparisons),
                   dimnames = list(row_label, selected_comparisons))
padj_mat <- lfc_mat

for (comp_name in selected_comparisons) {
  res_file <- file.path(deseq2_dir, paste0(comp_name, "_DESeq2_results.csv"))
  if (!file.exists(res_file)) stop("Missing DESeq2 results file: ", res_file)
  
  res_df <- read.csv(res_file, stringsAsFactors = FALSE)
  
  match_idx <- match(hsp_genes$GeneID, res_df$Gene)
  missing   <- is.na(match_idx)
  if (any(missing)) {
    fallback <- match(strip_version(hsp_genes$GeneID[missing]), strip_version(res_df$Gene))
    match_idx[missing] <- fallback
  }
  
  lfc_mat[, comp_name]  <- res_df$log2FoldChange[match_idx]
  padj_mat[, comp_name] <- res_df$padj[match_idx]
  
  n_found <- sum(!is.na(match_idx))
  message(comp_name, ": found ", n_found, " / ", nrow(hsp_genes), " HSP70 genes")
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
          file.path(output_dir, "HSP70_Heatmap_Column_Legend.csv"), row.names = FALSE)

# Swap in the short a-h labels for the plot itself
colnames(lfc_mat)  <- comparison_labels[colnames(lfc_mat)]
colnames(star_mat) <- comparison_labels[colnames(star_mat)]

italic_labels <- as.expression(lapply(rownames(lfc_mat), function(g) bquote(italic(.(g)))))

# symmetric diverging scale so 0 (no change) sits at white
max_abs <- max(abs(lfc_mat), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)
palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

tiff(file.path(output_dir, "HSP70_Heatmap.tiff"),   # saved directly into Results/7_HSP
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
  angle_col        = 0,           # single letters (a-h), horizontal reads fine and saves space
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

message("\nDone. Saved: ", normalizePath(file.path(output_dir, "HSP70_Heatmap.tiff")))