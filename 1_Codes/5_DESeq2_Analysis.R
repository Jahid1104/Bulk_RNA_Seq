# ==========================================================
# 5_DESeq2_Analysis.R by Md Jahid Hasan Jone
# Simple DESeq2 differential expression pipeline
# ==========================================================
# What this script does:
#   1. Reads the gene count csv and metadata csv (direct path, or pick them with file.choose(); see section 2)
#   2. Runs DESeq2 separately for each comparison listed in "comparisons" below
#   3. Saves a results table (csv) for every comparison
#   4. Saves Volcano plot, MA plot, and Heatmap for every comparison
#   5. Saves one combined csv of the top genes (by adjusted p-value) from every comparison
#   6. Saves one bar plot summarizing up/down regulated gene counts across all comparisons, 
#      and one Venn diagram (ggVennDiagram) covering every comparison's DEGs
#   7. Saves an overall sample PCA (with 95% ellipses) and a top-variable-genes heatmap (all samples)
#
# Sample_ID format reminder: 1F_T0_24h_R1
#   1/2 = Genotype | F/L = Tissue | T0/T1 = Temperature | 24h/48h/72h = Time | R# = Replicate
#
# ----------------------------------------------------------
# EXPERIMENTAL SETUP
# ----------------------------------------------------------
#   Genotype     1 = CLN1466EA, 2 = NC123S
#   Tissue       F = Flower, L = Leaf
#   Temperature  T0 = control, T1 = heat
#   Time         24h, 48h, 72h
#   Replicate    R1, R2, ... (biological replicates of each combination)
# Metadata.csv has one row per sample with the columns:
#   Sample_ID, Genotype, Tissue, Temperature, Time
# Temperature is not used in the comparison list below (the "Heat vs Control"
# line is commented out), so T0 and T1 samples are pooled inside each comparison.
#
# ----------------------------------------------------------
# HOW TO EDIT THE COMPARISON LIST (section 5)
# ----------------------------------------------------------
# Each comparison is one line:
#   list(name = "...", group_col = "...", group1 = "...", group2 = "...", filter = ...)
#     name       label used in file names, plots and tables; must be unique and
#                must not contain characters that are not allowed in file names (e.g. "/")
#     group_col  metadata column that holds the two groups: Genotype, Tissue, Temperature or Time
#     group1     value compared against group2 (positive log2 fold change = higher in group1)
#     group2     reference value
#     filter     NULL = use all samples, or keep only some samples, e.g.
#                list(Tissue = "F", Time = "24h") = Flower samples at 24h only.
#                Filter on any metadata column except group_col.
# Each comparison runs its own DESeq2 model (design = ~ group_col) on the samples that
# pass the filter. Factors that are in neither group_col nor filter are pooled.
#
# Add a comparison:    copy a line, change the fields, e.g.
#   list(name = "Heat_Flower vs Control_Flower", group_col = "Temperature", group1 = "T1", group2 = "T0", filter = list(Tissue = "F"))
# Remove a comparison: delete its line or put # in front of it.
# Every line needs a comma at the end except the last one in the list.
#
# After changing the list, also check:
#   - venn_comparisons (section 10) and final_comparisons (near the end) must use names
#     exactly as written in "name" above, or the Venn / combined volcano figures fail.
#   - The Venn diagram has 6 fill colors (set_colors), so it handles up to 6 comparisons.
#   - The supplementary MA and volcano figures are 5 x 5 grids (25 comparisons, letters a-y).
#     For more than 25 comparisons, increase ncol / nrow in the two wrap_plots() calls.
# ==========================================================





# -----------------------------
# 0. Packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
if (!requireNamespace("DESeq2",      quietly = TRUE)) 
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)

cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "ggVennDiagram", "extrafont", "patchwork")
to_install <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggVennDiagram)
library(extrafont)
library(patchwork)




# -----------------------------
# 1. Font (Times New Roman, gene names italic)
# -----------------------------
# First time EVER using extrafont on this computer, run this once (takes a few minutes), 
#then you can comment it out again:

#font_import()
loadfonts(device = if (.Platform$OS.type == "windows") "win" else "pdf", quiet = TRUE)   # Windows. Use device = "pdf" on Mac/Linux
FONT <- "Times New Roman"                 # <- EDIT: use "serif" here if this font isn't found

theme_set(theme_bw(base_family = FONT) + theme(text = element_text(family = FONT)))





# -----------------------------
# 2. Load data
# -----------------------------
# Input files  <- EDIT
# There are two ways to give each input file. Use ONE of them per file and comment out the other:
#   Option 1 (default): direct path. Replace /.../.../ with the path to your project folder.
#   Option 2: file.choose() opens a window to pick the file (needs an interactive R session such
#             as RStudio). To use it, remove the # from the "message" and "file.choose()" lines
#             and put a # in front of the direct-path line.


# Gene count csv (Gene_ID = 1st column, samples = other columns)
# Option 1: direct path
counts_file <- "/.../.../Bulk_RNA_seq/5_Results/3_tximport/gene_counts.csv"

# Option 2: choose the file in a window
#message("Choose the GENE COUNT csv file (Gene_ID = 1st column, samples = other columns)")
#counts_file <- file.choose()

# Metadata csv (Sample_ID = 1st column; columns Genotype, Tissue, Temperature, Time)
# Option 1: direct path
meta_file <- "/.../.../Bulk_RNA_seq/Metadata.csv"

# Option 2: choose the file in a window
#message("Choose the METADATA csv file (Sample_ID = 1st column)")
#meta_file <- file.choose()


counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(round(counts))          # DESeq2 needs integer counts
#View(counts)


meta <- read.csv(meta_file, row.names = "Sample_ID")
meta <- meta[colnames(counts), , drop = FALSE]   # align sample order to the count matrix
meta$Genotype <- factor(meta$Genotype)
#View(meta)


output_dir <- "/.../.../Bulk_RNA_seq/5_Results/5_DESeq2"   # <- EDIT: folder for all outputs
dir.create(file.path(output_dir, "Tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Tables/Comparisons"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Figures/Heatmaps"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Figures/Volcano_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Figures/MA_Plots"),      recursive = TRUE, showWarnings = FALSE)





# -----------------------------
# 3. Parameters  <- EDIT any of these
# -----------------------------
padj_cutoff   <- 0.05   # adjusted p-value (FDR) cutoff for "significant"
lfc_cutoff    <- 2      # log2 fold-change cutoff for "significant" (1 = 2-fold change)
min_gene_count<- 10     # drop genes with fewer total reads than this across the comparison's samples
top_n_labels  <- 10     # number of top genes labeled on each volcano plot
top_n_heatmap <- 50     # max number of genes shown per heatmap





# -----------------------------
# 4. Overall QC: PCA + heatmap of most variable genes
# -----------------------------
dds_all <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
dds_all <- dds_all[rowSums(counts(dds_all)) > 0, ]
vsd_all <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd_all, intgroup = c("Genotype", "Tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$Group <- paste0(pca_data$Genotype, pca_data$Tissue)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = name)) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 0.4) +  # 95% CI ellipse per group
  geom_point(size = 3) +
  geom_text_repel(size = 3, family = FONT, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle(" ")

ggsave(file.path(output_dir, "Figures", "PCA_AllSamples.tiff"), p_pca,
       dpi = 1000, width = 8, height = 6, compression = "lzw")

rv <- apply(assay(vsd_all), 1, var)
top_var_genes <- names(sort(rv, decreasing = TRUE))[1:top_n_heatmap]
mat_all <- assay(vsd_all)[top_var_genes, ]
mat_all <- mat_all - rowMeans(mat_all)
italic_all <- as.expression(lapply(rownames(mat_all), function(g) bquote(italic(.(g)))))


tiff(file.path(output_dir, "Figures/Heatmaps", "Heatmap_TopVariableGenes_AllSamples.tiff"),
     width = 10, height = 10, units = "in", res = 1000, compression = "lzw")

pheatmap(
  mat_all,
  annotation_col = meta[, c("Genotype", "Tissue", "Temperature", "Time")],
  labels_row = italic_all,
  fontsize_row = 6,
  fontsize_col = 6,       # <-- column/sample name font size
  fontfamily = FONT,
  show_colnames = TRUE
)

dev.off()





# -----------------------------
# 5. Define comparisons  <- EDIT to add/remove/change comparisons
# -----------------------------
# group_col  = metadata column that defines the two groups being compared
# group1     = "treatment" side of the comparison (positive log2FC = higher in group1)
# group2     = reference side
# filter     = optional extra restriction, e.g. list(Tissue = "F") to only use Flower samples
comparisons <- list(
  list(name = "CLN1466EA vs NC123S", group_col = "Genotype", group1 = "1", group2 = "2", filter = NULL),
  
  list(name = "Flower vs Leaf",      group_col = "Tissue",   group1 = "F", group2 = "L", filter = NULL),
  
  #list(name = "Heat vs Control",     group_col = "Temperature", group1 = "T1", group2 = "T0", filter = NULL),
  
  list(name = "24h vs 48h",          group_col = "Time",     group1 = "24h", group2 = "48h", filter = NULL),
  list(name = "24h vs 72h",          group_col = "Time",     group1 = "24h", group2 = "72h", filter = NULL),
  list(name = "48h vs 72h",          group_col = "Time",     group1 = "48h", group2 = "72h", filter = NULL),
  
  list(name = "CLN1466EA_Flower_24h vs CLN1466EA_Flower_48h", group_col = "Time", group1 = "24h", group2 = "48h", filter = list(Genotype = "1", Tissue = "F")),
  list(name = "CLN1466EA_Flower_24h vs CLN1466EA_Flower_72h", group_col = "Time", group1 = "24h", group2 = "72h", filter = list(Genotype = "1", Tissue = "F")),
  list(name = "CLN1466EA_Flower_48h vs CLN1466EA_Flower_72h", group_col = "Time", group1 = "48h", group2 = "72h", filter = list(Genotype = "1", Tissue = "F")),
  
  list(name = "CLN1466EA_Leaf_24h vs CLN1466EA_Leaf_48h", group_col = "Time", group1 = "24h", group2 = "48h", filter = list(Genotype = "1", Tissue = "L")),
  list(name = "CLN1466EA_Leaf_24h vs CLN1466EA_Leaf_72h", group_col = "Time", group1 = "24h", group2 = "72h", filter = list(Genotype = "1", Tissue = "L")),
  list(name = "CLN1466EA_Leaf_48h vs CLN1466EA_Leaf_72h", group_col = "Time", group1 = "48h", group2 = "72h", filter = list(Genotype = "1", Tissue = "L")),
  
  list(name = "NC123S_Flower_24h vs NC123S_Flower_48h", group_col = "Time", group1 = "24h", group2 = "48h", filter = list(Genotype = "2", Tissue = "F")),
  list(name = "NC123S_Flower_24h vs NC123S_Flower_72h", group_col = "Time", group1 = "24h", group2 = "72h", filter = list(Genotype = "2", Tissue = "F")),
  list(name = "NC123S_Flower_48h vs NC123S_Flower_72h", group_col = "Time", group1 = "48h", group2 = "72h", filter = list(Genotype = "2", Tissue = "F")),
  
  list(name = "NC123S_Leaf_24h vs NC123S_Leaf_48h", group_col = "Time", group1 = "24h", group2 = "48h", filter = list(Genotype = "2", Tissue = "L")),
  list(name = "NC123S_Leaf_24h vs NC123S_Leaf_72h", group_col = "Time", group1 = "24h", group2 = "72h", filter = list(Genotype = "2", Tissue = "L")),
  list(name = "NC123S_Leaf_48h vs NC123S_Leaf_72h", group_col = "Time", group1 = "48h", group2 = "72h", filter = list(Genotype = "2", Tissue = "L")), 
  
  list(name = "CLN1466EA_Flower vs NC123S_Flower", group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "F")),
  
  list(name = "CLN1466EA_Flower_24h vs NC123S_Flower_24h", group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "F", Time = "24h")),
  list(name = "CLN1466EA_Flower_48h vs NC123S_Flower_48h", group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "F", Time = "48h")),
  list(name = "CLN1466EA_Flower_72h vs NC123S_Flower_72h", group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "F", Time = "72h")),
  
  list(name = "CLN1466EA_Leaf vs NC123S_Leaf",      group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "L")),
  
  list(name = "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h",     group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "L", Time = "24h")),
  list(name = "CLN1466EA_Leaf_48h vs NC123S_Leaf_48h",     group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "L", Time = "48h")),
  list(name = "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h",     group_col = "Genotype", group1 = "1", group2 = "2", filter = list(Tissue = "L", Time = "72h"))
)





# -----------------------------
# 6. Functions (run once per comparison, don't need to edit)
# -----------------------------
run_comparison <- function(comp) {
  sub_meta <- meta
  if (!is.null(comp$filter)) {
    for (col in names(comp$filter)) sub_meta <- sub_meta[sub_meta[[col]] == comp$filter[[col]], ]
  }
  sub_meta <- sub_meta[sub_meta[[comp$group_col]] %in% c(comp$group1, comp$group2), , drop = FALSE]
  sub_meta[[comp$group_col]] <- factor(sub_meta[[comp$group_col]], levels = c(comp$group2, comp$group1))
  sub_counts <- counts[, rownames(sub_meta)]
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_meta,
                                design = as.formula(paste0("~", comp$group_col)))
  dds <- dds[rowSums(counts(dds)) >= min_gene_count, ]
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c(comp$group_col, comp$group1, comp$group2), alpha = padj_cutoff)
  res <- lfcShrink(dds, contrast = c(comp$group_col, comp$group1, comp$group2), res = res, type = "normal")
  
  res_df <- as.data.frame(res)
  res_df$Gene <- rownames(res_df)
  res_df <- res_df[, c("Gene", setdiff(colnames(res_df), "Gene"))]
  res_df$Regulation <- "NS"
  res_df$Regulation[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff] <- "Up"
  res_df$Regulation[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down"
  
  write.csv(res_df, file.path(output_dir, "Tables/Comparisons", paste0(comp$name, "_DESeq2_results.csv")), row.names = FALSE)
  list(dds = dds, res_df = res_df, name = comp$name, group_col = comp$group_col)
}



make_volcano <- function(res_df, comp_name) {
  res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "NS"))
  top_genes <- res_df[!is.na(res_df$padj) & res_df$Regulation != "NS", ]
  top_genes <- head(top_genes[order(top_genes$padj), ], top_n_labels)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, fontface = "italic",
                    max.overlaps = 20, show.legend = FALSE) +
    labs(title = " ", x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  
  ggsave(file.path(output_dir, "Figures/Volcano_Plots", paste0("Volcano_", comp_name, ".tiff")),
         p, dpi = 1000, width = 8, height = 8, compression = "lzw")
}



make_ma_plot <- function(res_df, comp_name) {
  res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "NS"))
  
  p <- ggplot(res_df, aes(x = log2(baseMean + 1), 
                          y = log2FoldChange, 
                          color = Regulation)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = " ", x = "log2(mean expression + 1)", y = "log2 Fold Change")
  
  ggsave(file.path(output_dir, "Figures/MA_Plots", paste0("MAplot_", comp_name, ".tiff")),
         p, dpi = 1000, width = 8, height = 8, compression = "lzw")
}




make_comparison_heatmap <- function(dds, res_df, comp_name, group_col) {
  sig_genes <- res_df$Gene[res_df$Regulation != "NS"]
  if (length(sig_genes) < 2) { message("  (skipping heatmap for ", comp_name, ": fewer than 2 DEGs)"); return(invisible(NULL)) }
  
  vsd <- vst(dds, blind = FALSE)
  ordered_sig <- res_df$Gene[order(res_df$padj)]
  ordered_sig <- ordered_sig[ordered_sig %in% sig_genes]
  keep <- head(ordered_sig, top_n_heatmap)
  
  mat <- assay(vsd)[keep, ]
  mat <- mat - rowMeans(mat)
  italic_labels <- as.expression(lapply(rownames(mat), function(g) bquote(italic(.(g)))))
  annotation_col <- as.data.frame(colData(dds)[, group_col, drop = FALSE])
  
  tiff(file.path(output_dir, "Figures/Heatmaps", paste0("Heatmap_", comp_name, ".tiff")),
       width = 10, height = 10, units = "in", res = 1000, compression = "lzw")
  
  pheatmap(mat, 
           annotation_col = annotation_col, 
           labels_row = italic_labels,
           fontsize_row = 7, 
           fontsize_col = 6, 
           fontfamily = FONT, 
           show_colnames = TRUE)
  
  dev.off()
}





# -----------------------------
# 7. Run every comparison
# -----------------------------
all_results <- list()
n_comparisons <- length(comparisons)

for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  message(sprintf("[%d/%d] Running comparison: %s ...", i, n_comparisons, comp$name))
  out <- run_comparison(comp)
  make_volcano(out$res_df, out$name)
  make_ma_plot(out$res_df, out$name)
  make_comparison_heatmap(out$dds, out$res_df, out$name, out$group_col)
  all_results[[comp$name]] <- out$res_df
  message(sprintf("[%d/%d] Completed: %s", i, n_comparisons, comp$name))
}

message("\nAll ", n_comparisons, " comparisons completed.")


# -----------------------------
# 8a. Top genes per comparison (combined csv)
# -----------------------------
top_n_export <- 10   # <- EDIT: how many top genes (by lowest adjusted p-value) to export per comparison
top_genes_all <- do.call(rbind, lapply(names(all_results), function(nm) {
  df <- all_results[[nm]]
  df <- df[!is.na(df$padj) & df$Regulation != "NS", ]
  if (nrow(df) == 0) return(NULL)   # skip comparisons with no significant genes
  df <- df[order(df$padj), ]
  df <- head(df, top_n_export)
  cbind(Comparison = nm, df)
}))
write.csv(top_genes_all,
          file.path(output_dir, "Tables", paste0("Top", top_n_export, "_Genes_PerComparison.csv")),
          row.names = FALSE)

# -----------------------------
# 8b. Top 10 DE genes overall (for qRT-PCR validation)
# -----------------------------
# Ranks genes by how many comparisons they are significant DE in, then by
# average |log2FoldChange| across those comparisons -- genes that are
# consistently and strongly DE make the best qPCR validation candidates.

gene_summary <- do.call(rbind, lapply(names(all_results), function(nm) {
  df <- all_results[[nm]]
  df <- df[df$Regulation != "NS", c("Gene", "log2FoldChange", "padj")]
  if (nrow(df) == 0) return(NULL)
  df$Comparison <- nm
  df
}))

top_de_overall <- do.call(rbind, lapply(split(gene_summary, gene_summary$Gene), function(g) {
  data.frame(
    Gene = g$Gene[1],
    N_Comparisons_Significant = nrow(g),
    Mean_abs_log2FC = mean(abs(g$log2FoldChange)),
    Min_padj = min(g$padj),
    Comparisons = paste(g$Comparison, collapse = "; ")
  )
}))

top_de_overall <- top_de_overall[order(-top_de_overall$N_Comparisons_Significant, -top_de_overall$Mean_abs_log2FC), ]
top_de_overall <- head(top_de_overall, 10)   # <- EDIT: change 10 to however many candidates you want

write.csv(top_de_overall, file.path(output_dir, "Tables", "Top10_DE_Genes_Overall_for_qPCR.csv"), row.names = FALSE)





# -----------------------------
# 9. Bar plot: up/down regulated genes per comparison
# -----------------------------

summary_df <- do.call(rbind, lapply(names(all_results), function(nm) {
  df <- all_results[[nm]]
  data.frame(Comparison = nm, Up = sum(df$Regulation == "Up"), Down = sum(df$Regulation == "Down"))
}))
write.csv(summary_df, file.path(output_dir, "Tables", "DEG_Summary_UpDown.csv"), row.names = FALSE)

summary_long <- rbind(
  data.frame(Comparison = summary_df$Comparison, Regulation = "Up",   Count =  summary_df$Up),
  data.frame(Comparison = summary_df$Comparison, Regulation = "Down", Count = -summary_df$Down)
)
summary_long$Comparison <- factor(summary_long$Comparison, levels = summary_df$Comparison)  # keep original comparison order

p_bar <- ggplot(summary_long, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_col() +
  geom_text(data = subset(summary_long, Count != 0),   # skip "0" labels -- no bar to attach to, and they collide with the other label
            aes(label = abs(Count), vjust = ifelse(Count >= 0, -0.4, 1.4)),
            size = 3, family = FONT, show.legend = FALSE) +   # <- EDIT: change "size" to make bar labels bigger/smaller
  scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  geom_hline(yintercept = 0) +
  labs(title = "Number of up- and down-regulated genes", x = "", y = "Number of genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(output_dir, "Figures", "DEG_Summary_BarPlot.tiff"),
       p_bar, dpi = 1000, width = 10, height = 8, compression = "lzw")





# -----------------------------
# 10. Venn diagram (six selected comparisons, single diagram)
# -----------------------------
# https://gaospecial.github.io/ggVennDiagram/
# ggVennDiagram natively draws 2-7 sets as circles/ellipses/polygons. With
# more than 7 sets it switches to an upset plot instead, which is
# a composite of several sub-plots rather than one ggplot object -- so the
# fill/theme layers below only get added when a real Venn is being drawn.

venn_comparisons <- c(
  "CLN1466EA_Flower_24h vs NC123S_Flower_24h",
  "CLN1466EA_Flower_72h vs NC123S_Flower_72h",
  "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h",
  "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h",
  "CLN1466EA_Flower vs NC123S_Flower",   # genotype effect, Flower only, pooled across time
  "CLN1466EA_Leaf vs NC123S_Leaf"        # genotype effect, Leaf only, pooled across time
)
venn_gene_lists <- lapply(all_results[venn_comparisons], function(df) {
  df$Gene[!is.na(df$padj) & df$Regulation != "NS"]
})

# Build the Venn geometry manually so each comparison's circle gets its own flat
# background color instead of a count-based gradient fill
venn_obj  <- Venn(venn_gene_lists)
venn_data <- process_data(venn_obj)

set_edges   <- venn_setedge(venn_data)      # circle/ellipse outline per comparison (columns: X, Y, id)
set_labels  <- venn_setlabel(venn_data)     # label position per comparison (columns: X, Y, id, name)
region_lbls <- venn_regionlabel(venn_data)  # intersection count label positions (columns: X, Y, count)

id_to_name <- unique(set_labels[, c("id", "name")])
set_edges  <- merge(set_edges, id_to_name, by = "id", sort = FALSE)

# a, b, c... shown on the diagram; full names shown only in the legend
short_codes <- setNames(LETTERS[seq_along(venn_comparisons)], venn_comparisons)
set_colors  <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")[seq_along(venn_comparisons)],  # <- EDIT: change hex colors here
  venn_comparisons
)
set_labels$short   <- short_codes[set_labels$name]
legend_labels       <- paste0(short_codes, ": ", names(short_codes))
names(legend_labels) <- names(short_codes)

p_venn <- ggplot() +
  geom_polygon(data = set_edges, aes(X, Y, group = id, fill = name), color = "black", alpha = 0.4) +
  geom_label(data = region_lbls, aes(X, Y, label = count), size = 3, family = FONT, label.size = 0) +
  geom_text(data = set_labels, aes(X, Y, label = short), size = 5, fontface = "bold", family = FONT) +
  scale_fill_manual(
    values = set_colors,
    breaks = names(set_colors),
    labels = legend_labels[names(set_colors)],
    name = "Comparison"
  ) +
  coord_equal() +
  theme_void(base_family = FONT) +
  theme(legend.position = "right", legend.text = element_text(size = 8, family = FONT))

ggsave(file.path(output_dir, "Figures", "Venn_6Comparisons.tiff"),
       p_venn, dpi = 1000, width = 14, height = 10, compression = "lzw")





# -----------------------------
# Supplementary Figure: combined MA plots (5x5 grid, A4-sized)
# -----------------------------

ma_letters <- letters[seq_along(all_results)]   # a, b, c, ... y (25 comparisons -> a to y)

make_ma_plot_mini <- function(res_df, panel_letter) {
  res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "NS"))
  ggplot(res_df, aes(x = log2(baseMean + 1), y = log2FoldChange, color = Regulation)) +
    geom_point(alpha = 0.6, size = 0.4) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(title = panel_letter, x = NULL, y = NULL) +
    theme_bw(base_family = FONT, base_size = 8) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 6))
}

ma_plot_list <- lapply(seq_along(all_results), function(i) make_ma_plot_mini(all_results[[i]], ma_letters[i]))

combined_ma <- wrap_plots(ma_plot_list, ncol = 5, nrow = 5) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(family = FONT, size = 14, face = "bold", hjust = 0.5))
  )

# Panel letter -> full comparison name, so you know which letter is which when you add labels
ma_panel_legend <- data.frame(Panel = ma_letters, Comparison = names(all_results))
write.csv(ma_panel_legend, file.path(output_dir, "Tables", "MA_Panel_Letter_Legend.csv"), row.names = FALSE)

# Sized to A4 (8.27 x 11.69 in) minus a 1-inch margin on every side, portrait
# orientation -> fits a printed A4 page at 100% scale with no resizing needed.
# For landscape instead, swap width/height (9.69 x 6.27).
ggsave(file.path(output_dir, "Figures", "Supplementary_MA_AllComparisons.tiff"),
       combined_ma, dpi = 1000, width = 6.27, height = 9.69, units = "in", compression = "lzw")


# -----------------------------
# Supplementary Figure: combined volcano plots (5x5 grid, A4-sized)
# -----------------------------
volcano_letters <- letters[seq_along(all_results)]   # a, b, c, ... y (25 comparisons -> a to y)

make_volcano_mini <- function(res_df, panel_letter) {
  res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "NS"))
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.6, size = 0.4) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.3) +
    labs(title = panel_letter, x = NULL, y = NULL) +
    theme_bw(base_family = FONT, base_size = 8) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.position = "none",
          axis.text = element_text(size = 6))
}

volcano_plot_list <- lapply(seq_along(all_results), function(i) make_volcano_mini(all_results[[i]], volcano_letters[i]))

combined_volcano <- wrap_plots(volcano_plot_list, ncol = 5, nrow = 5) +
  plot_annotation(
    title = " ",
    theme = theme(plot.title = element_text(family = FONT, size = 14, face = "bold", hjust = 0.5))
  )

# Panel letter -> full comparison name
volcano_panel_legend <- data.frame(Panel = volcano_letters, Comparison = names(all_results))
write.csv(volcano_panel_legend, file.path(output_dir, "Tables", "Volcano_Panel_Letter_Legend.csv"), row.names = FALSE)

# A4 portrait minus 1-inch margins, same as the MA plot grid
ggsave(file.path(output_dir, "Figures", "Supplementary_Volcano_AllComparisons.tiff"),
       combined_volcano, dpi = 1000, width = 6.27, height = 9.69, units = "in", compression = "lzw")


# -----------------------------
# Combined Volcano Figure (panels a-h volcano)
# -----------------------------
# Layout: 2 rows x 4 columns
#   Row 1: a  b  c  d
#   Row 2: e  f  g  h

final_comparisons <- c(
  "CLN1466EA_Flower_24h vs NC123S_Flower_24h",
  "CLN1466EA_Flower_72h vs NC123S_Flower_72h",
  "CLN1466EA_Leaf_24h vs NC123S_Leaf_24h",
  "CLN1466EA_Leaf_72h vs NC123S_Leaf_72h",
  "CLN1466EA_Flower vs NC123S_Flower",
  "CLN1466EA_Leaf vs NC123S_Leaf",
  "CLN1466EA vs NC123S",
  "Flower vs Leaf"
)
top_n_labels_final <- 5   # <- EDIT: number of top genes labeled on each panel a-h

make_volcano_final <- function(res_df) {
  res_df$Regulation <- factor(res_df$Regulation, levels = c("Up", "Down", "NS"))
  top_genes <- res_df[!is.na(res_df$padj) & res_df$Regulation != "NS", ]
  top_genes <- head(top_genes[order(top_genes$padj), ], top_n_labels_final)
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Regulation)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.3) +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 1.8, fontface = "italic",
                    family = FONT, max.overlaps = 20, segment.linewidth = 0.2, show.legend = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_family = FONT, base_size = 8) +
    theme(legend.position = "none", axis.text = element_text(size = 4))
}

final_volcano_list <- lapply(final_comparisons, function(nm) make_volcano_final(all_results[[nm]]))

combined_final <- wrap_plots(final_volcano_list, ncol = 4, nrow = 2) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(plot.tag = element_text(family = FONT, size = 14, face = "bold"))
  )

# Panel letter -> comparison name, for reference
final_panel_legend <- data.frame(Panel = letters[1:8], Comparison = final_comparisons)
write.csv(final_panel_legend, file.path(output_dir, "Tables", "Combined_Volcano_Figure_Panel_Legend.csv"), row.names = FALSE)

ggsave(file.path(output_dir, "Figures", "Combined_Volcano_Figure.tiff"),
       combined_final, dpi = 1000, width = 6.27, height = 8, units = "in", compression = "lzw")


message("\nDone. All tables and figures saved to: ", normalizePath(output_dir))