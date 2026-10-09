# ==========================================================
# 8_WGCNA.R by Md Jahid Hasan Jone
# Weighted gene co-expression network analysis (WGCNA)
# ==========================================================
# What this script does:
#   Reads the gene count matrix from step 3.2 and Metadata.csv, removes low-expression genes,
#   applies a variance stabilizing transformation (VST), keeps the most variable genes and builds
#   a signed co-expression network. Every panel is saved standalone (1000 dpi TIFF):
#
#   A.  Soft threshold (scale independence + mean connectivity)
#   B.  Gene dendrogram, modules, + per-treatment-group correlation rows
#   C.  Eigengene dendrogram + eigengene adjacency heatmap
#   D.  TOM network heatmap (random gene subsample)
#   E.  Module-trait heatmap, factor-level traits (Genotype/Tissue/Temp/Time)
#   F.  Module-group heatmap, one column per exact treatment group
#   G.  Module membership vs. gene significance (Genotype within Flower/Leaf)
#   H.  Module expression heatmap + eigengene barplot, per focal module
#   I.  Module co-expression network graph, hub genes in red
#
# Focal modules for panels G and H are the modules most correlated with
# Genotype-within-Flower and Genotype-within-Leaf. Panel I uses the modules
# listed in "panel_I_modules" (section 23).
#
# Table 1 equivalent: Hub_Gene_Annotation_Table.csv. It needs a gene annotation file
# (the GFF, see section 3). Table 1 always lists exactly the hub genes shown in Panel I
# (same modules, same top-5 genes).
#
# Panel I combines its module figures into one 2x2 image with the combine_panels_grid()
# helper in section 3 (see section 23.3b). Text is set to a minimum of 10 pt throughout.
#
# Metadata.csv must have Sample_ID as the 1st column and the columns Genotype, Tissue,
# Temperature and Time (Tissue values F = Flower, L = Leaf).
#
# ----------------------------------------------------------
# HOW TO EDIT THE SCRIPT
# ----------------------------------------------------------
#   Paths (section 3): replace /.../.../ with the path to your project folder. Each input
#     file can also be picked in a window with file.choose() (see the Option 2 lines).
#   Network settings (section 3): number of genes, filtering, module size, merge height, network type.
#   Modules in Panel I (section 23): module colors depend on your data. Run the script once,
#     check Module_Sizes.csv, then edit "panel_I_modules". The 2x2 figure in section 23.3b
#     is written for exactly 4 modules; edit that part too if you change the number.
#   Different experiment: the traits in section 10 use the metadata columns Genotype, Tissue,
#     Temperature and Time (and Tissue = F / L). Edit section 10 if your metadata differs.
# ==========================================================


# -----------------------------
# 1. Clean environment
# -----------------------------

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)


# -----------------------------
# 2. Packages
# -----------------------------

cran_packages <- c("WGCNA", "ggplot2", "pheatmap", "RColorBrewer", "igraph", "magick")
bioc_packages <- c("DESeq2", "impute", "preprocessCore")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

library(WGCNA)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(magick)

enableWGCNAThreads()
options(stringsAsFactors = FALSE)

# DESeq2's Bioconductor dependencies redefine cor() after WGCNA is loaded,
# which drops arguments (weights.x, weights.y, cosine) that blockwiseModules,
# pickSoftThreshold, and adjacency() expect. Force cor() back to WGCNA's version.
cor <- WGCNA::cor


# -----------------------------
# 3. User settings
# -----------------------------

# Input files  <- EDIT
# There are two ways to give each input file. Use ONE of them per file and comment out the other:
#   Option 1 (default): direct path. Replace /.../.../ with the path to your project folder.
#   Option 2: file.choose() opens a window to pick the file. To use it, remove the # from the
#             "message" and "file.choose()" lines and put a # in front of the direct-path line.

# Gene count csv from step 3.2 (Gene_ID = 1st column, samples = other columns)
# Option 1: direct path
count_file    <- "/.../.../Bulk_RNA_seq/5_Results/3_tximport/gene_counts.csv"

# Option 2: choose the file in a window
#message("Choose the GENE COUNT csv file (Gene_ID = 1st column, samples = other columns)")
#count_file <- file.choose()

# Metadata csv (Sample_ID = 1st column; columns Genotype, Tissue, Temperature, Time)
# Option 1: direct path
metadata_file <- "/.../.../Bulk_RNA_seq/Metadata.csv"

# Option 2: choose the file in a window
#message("Choose the METADATA csv file (Sample_ID = 1st column)")
#metadata_file <- file.choose()

# GFF3 file (ITAG4.0), used to build the Table 1 hub gene annotation from the mRNA "Note=" field.
# If the file is not found, the rest of the analysis still runs and only Table 1 is skipped.
# Option 1: direct path
gff_file <- "/.../.../Bulk_RNA_seq/2_References/gene_annotation.gff"

# Option 2: choose the file in a window
#message("Choose the GFF3 annotation file (gene_annotation.gff)")
#gff_file <- file.choose()

# Output folder  <- EDIT
output_dir  <- "/.../.../Bulk_RNA_seq/5_Results/8_WGCNA"
figures_dir <- file.path(output_dir, "Figures")

if (!dir.exists(output_dir))  dir.create(output_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# Network settings  <- EDIT if needed
N_GENES          <- 5000      # number of most variable genes used to build the network (a smaller value is used if fewer pass the filter)
MIN_CPM          <- 1         # a gene is kept if it has at least this CPM ...
MIN_SAMPLE_PROP  <- 0.20      # ... in at least this fraction of samples (and in at least 3 samples)
MIN_MODULE_SIZE  <- 30        # smallest module, in genes
MERGE_CUT_HEIGHT <- 0.25      # modules whose eigengenes are closer than this are merged
NETWORK_TYPE     <- "signed"  # "signed" or "unsigned"
CORRELATION_FUNCTION <- "bicor"   # "bicor" or "pearson"
KME_THRESHOLD    <- 0.80      # a gene is saved as a hub gene if its absolute module membership |kME| is at least this

# Plot settings
N_GENES_TOM_PLOT <- 1000      # random gene subsample drawn in the TOM heatmap (panel D)
EDGE_WEIGHT_THRESHOLD <- 0.15 # weakest connection drawn in the network graphs (panel I)

FIG_DPI       <- 1000
FIG_POINTSIZE <- 10
FIG_COMPRESSION <- "lzw"

save_panel_tiff <- function(filename, width_in, height_in, plot_fun) {
  tiff(filename, width = width_in, height = height_in, units = "in",
       res = FIG_DPI, pointsize = FIG_POINTSIZE, compression = FIG_COMPRESSION)
  plot_fun()
  dev.off()
}

# combine_panels_to_figure() -- assembles a multi-panel figure by stacking
# ALREADY-RENDERED, standalone panel TIFFs (read from disk with magick and
# resized to a common width), rather than re-calling the panel plot
# functions inside a shared layout()/par(mfrow) block.
#
# Several panel functions (plotDendroAndColors for B, TOMplot for D) call
# layout()/par() internally, and panel F needs wide margins for its many
# treatment-group labels. Inside one shared layout these calls overwrite each
# other or fail with "figure margins too large". Each panel is therefore drawn
# to its own device first and the finished images are combined afterwards.
combine_panels_to_figure <- function(panel_files, out_file, target_width_px = 3000) {
  missing_files <- panel_files[!file.exists(panel_files)]
  if (length(missing_files) > 0) {
    stop("combine_panels_to_figure(): missing panel file(s): ", paste(missing_files, collapse = ", "))
  }
  imgs <- lapply(panel_files, image_read)
  imgs <- lapply(imgs, image_resize, geometry = paste0(target_width_px, "x"))
  combined <- image_append(image_join(imgs), stack = TRUE)
  image_write(combined, path = out_file, format = "tiff",
              density = paste0(FIG_DPI, "x", FIG_DPI))
  invisible(combined)
}

# combine_panels_grid() -- same idea as combine_panels_to_figure(), but for
# a 2D grid instead of a single stack. panel_files_matrix is a list of rows,
# each row a character vector of file paths (left to right). Every image in
# a row is resized to a common height first (so the row lines up), then rows
# are stacked and resized to a common width (so the columns line up). Used
# to combine Panel I's 4 module figures into one 2x2 image.
combine_panels_grid <- function(panel_files_matrix, out_file,
                                target_row_height_px = 1200) {
  all_files <- unlist(panel_files_matrix)
  missing_files <- all_files[!file.exists(all_files)]
  if (length(missing_files) > 0) {
    stop("combine_panels_grid(): missing panel file(s): ", paste(missing_files, collapse = ", "))
  }
  row_images <- lapply(panel_files_matrix, function(row_files) {
    imgs <- lapply(row_files, image_read)
    imgs <- lapply(imgs, image_resize, geometry = paste0("x", target_row_height_px))
    image_append(image_join(imgs), stack = FALSE)
  })
  grid <- image_append(image_join(row_images), stack = TRUE)
  image_write(grid, path = out_file, format = "tiff",
              density = paste0(FIG_DPI, "x", FIG_DPI))
  invisible(grid)
}


# -----------------------------
# 4. Read raw count matrix
# -----------------------------

counts <- read.csv(count_file, header = TRUE, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
mode(counts) <- "numeric"

if (anyNA(counts))   stop("NA values detected in the count matrix.")
if (any(counts < 0)) stop("Negative count values detected.")
if (anyDuplicated(rownames(counts)) > 0) counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]


# -----------------------------
# 5. Read metadata and match samples
# -----------------------------

metadata <- read.csv(metadata_file, header = TRUE, row.names = 1, check.names = FALSE)
metadata$Genotype <- factor(metadata$Genotype)
count_samples <- colnames(counts)

if (length(setdiff(count_samples, rownames(metadata))) > 0) stop("Sample mismatch between counts and metadata.")
if (length(setdiff(rownames(metadata), count_samples)) > 0) stop("Sample mismatch between counts and metadata.")

metadata <- metadata[count_samples, , drop = FALSE]
if (!identical(rownames(metadata), colnames(counts))) stop("Metadata sample order does not match count matrix.")

write.csv(metadata, file.path(output_dir, "Metadata_used_for_WGCNA.csv"))


# -----------------------------
# 6. Filter low-expression genes
# -----------------------------

library_sizes <- colSums(counts)
cpm_matrix <- sweep(counts, 2, library_sizes, "/") * 1e6
min_samples <- max(3, ceiling(ncol(counts) * MIN_SAMPLE_PROP))

keep_genes <- rowSums(cpm_matrix >= MIN_CPM) >= min_samples
counts_filtered <- counts[keep_genes, , drop = FALSE]
write.csv(counts_filtered, file.path(output_dir, "Filtered_Counts.csv"))


# -----------------------------
# 7. DESeq2 + VST (with local-fit fallback for large multi-factor datasets)
# -----------------------------

for (i in seq_len(ncol(metadata))) {
  if (is.character(metadata[[i]])) metadata[[i]] <- factor(metadata[[i]])
}

counts_filtered <- round(counts_filtered)
storage.mode(counts_filtered) <- "integer"

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata, design = ~ 1)

vsd <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) {
    message("vst() failed, falling back to local-fit VST: ", conditionMessage(e))
    dds_local <- estimateSizeFactors(dds)
    dds_local <- estimateDispersions(dds_local, fitType = "local")
    varianceStabilizingTransformation(dds_local, blind = TRUE)
  }
)

vst_matrix <- assay(vsd)
write.csv(vst_matrix, file.path(output_dir, "VST_Expression_Matrix.csv"))


# -----------------------------
# 8. Select most variable genes
# -----------------------------

gene_variance <- apply(vst_matrix, 1, var, na.rm = TRUE)
variance_table <- data.frame(Gene = rownames(vst_matrix), Variance = gene_variance)
variance_table <- variance_table[order(variance_table$Variance, decreasing = TRUE), ]

N_GENES <- min(N_GENES, nrow(variance_table))
selected_genes <- variance_table$Gene[seq_len(N_GENES)]
expression_wgcna <- vst_matrix[selected_genes, , drop = FALSE]

write.csv(variance_table[match(selected_genes, variance_table$Gene), ],
          file.path(output_dir, "WGCNA_Selected_Genes.csv"), row.names = FALSE)


# -----------------------------
# 9. Build WGCNA expression matrix and QC
# -----------------------------

datExpr <- as.data.frame(t(expression_wgcna))
gsg <- goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  datExpr  <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  metadata <- metadata[rownames(datExpr), , drop = FALSE]
}

pdf(file.path(output_dir, "01_Sample_Clustering.pdf"), width = 12, height = 8)
par(mar = c(5, 4, 4, 2))
plot(hclust(dist(datExpr), method = "average"), main = "Sample Clustering", xlab = "", sub = "", cex = 0.8)
dev.off()

pdf(file.path(output_dir, "02_Sample_Correlation_Heatmap.pdf"), width = 10, height = 9)
pheatmap(cor(t(datExpr), method = "pearson"), main = "Sample Correlation",
         fontsize = 10, fontsize_row = 10, fontsize_col = 10)
dev.off()


# -----------------------------
# 10. Build trait matrices
#     trait_matrix          -- one column per factor level + the two
#                               Genotype-within-tissue traits
#     trait_matrix_combined -- one column per exact treatment group
#                               (Genotype x Tissue x Temperature x Time),
#                               used for panels B and F (for example 1_F_T0_24h)
# -----------------------------

trait_matrix_list <- list()
for (trait_name in colnames(metadata)) {
  variable <- metadata[, trait_name]
  if (is.factor(variable) || is.character(variable)) {
    variable <- factor(variable)
    dummy_matrix <- model.matrix(~ variable - 1)
    colnames(dummy_matrix) <- paste0(trait_name, "_", levels(variable))
    rownames(dummy_matrix) <- rownames(metadata)
    trait_matrix_list[[trait_name]] <- dummy_matrix
  } else {
    trait_matrix_list[[trait_name]] <- matrix(as.numeric(variable), ncol = 1,
                                              dimnames = list(rownames(metadata), trait_name))
  }
}

trait_matrix <- as.data.frame(do.call(cbind, trait_matrix_list))
trait_matrix <- trait_matrix[rownames(datExpr), , drop = FALSE]

stopifnot(all(c("Genotype", "Tissue") %in% colnames(metadata)))

genotype_num <- setNames(as.numeric(factor(metadata$Genotype)) - 1, rownames(metadata))

genotype_in_flower <- setNames(rep(NA_real_, nrow(metadata)), rownames(metadata))
genotype_in_flower[metadata$Tissue == "F"] <- genotype_num[metadata$Tissue == "F"]

genotype_in_leaf <- setNames(rep(NA_real_, nrow(metadata)), rownames(metadata))
genotype_in_leaf[metadata$Tissue == "L"] <- genotype_num[metadata$Tissue == "L"]

trait_matrix$Genotype_within_Flower <- genotype_in_flower[rownames(trait_matrix)]
trait_matrix$Genotype_within_Leaf   <- genotype_in_leaf[rownames(trait_matrix)]

write.csv(trait_matrix, file.path(output_dir, "WGCNA_Trait_Matrix.csv"))

required_traits <- c("Genotype", "Tissue", "Temperature", "Time")
stopifnot(all(required_traits %in% colnames(metadata)))

metadata$Combined_Treatment <- paste(metadata$Genotype, metadata$Tissue, metadata$Temperature, metadata$Time, sep = "_")
combined_dummy <- model.matrix(~ Combined_Treatment - 1, data = metadata)
colnames(combined_dummy) <- gsub("^Combined_Treatment_?", "", colnames(combined_dummy))
rownames(combined_dummy) <- rownames(metadata)
trait_matrix_combined <- as.data.frame(combined_dummy[rownames(datExpr), , drop = FALSE])
write.csv(trait_matrix_combined, file.path(output_dir, "Combined_Experimental_Trait_Matrix.csv"))

# Sample order used by panels B, F, H: grouped by Genotype > Tissue > Temperature > Time
sample_order <- order(metadata$Genotype, metadata$Tissue, metadata$Temperature, metadata$Time)
group_order  <- unique(metadata$Combined_Treatment[sample_order])


# -----------------------------
# 11. Panel A - soft-thresholding power
# -----------------------------

powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = NETWORK_TYPE,
                         corFnc = CORRELATION_FUNCTION, verbose = 5)
write.csv(sft$fitIndices, file.path(output_dir, "Soft_Thresholding_Results.csv"), row.names = FALSE)

fit_indices <- sft$fitIndices
candidate_powers <- fit_indices$Power[fit_indices$SFT.R.sq >= 0.80]
softPower <- if (length(candidate_powers) > 0) candidate_powers[1] else fit_indices$Power[which.max(fit_indices$SFT.R.sq)]

plot_panel_A1 <- function() {
  plot(fit_indices[, 1], -sign(fit_indices[, 3]) * fit_indices[, 2], type = "n",
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology R^2", main = "Scale Independence")
  text(fit_indices[, 1], -sign(fit_indices[, 3]) * fit_indices[, 2], labels = fit_indices[, 1], cex = 1, col = "red")
  abline(h = 0.80, col = "blue", lty = 2)
}
plot_panel_A2 <- function() {
  plot(fit_indices[, 1], fit_indices[, 5], type = "n",
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", main = "Mean Connectivity")
  text(fit_indices[, 1], fit_indices[, 5], labels = fit_indices[, 1], cex = 1, col = "red")
}
save_panel_tiff(file.path(figures_dir, "Panel_A_Soft_Thresholding.tiff"), 10, 5,
                function() { par(mfrow = c(1, 2)); plot_panel_A1(); plot_panel_A2() })


# -----------------------------
# 12. Network construction and modules
# -----------------------------

net <- blockwiseModules(datExpr, power = softPower, TOMType = NETWORK_TYPE,
                        minModuleSize = MIN_MODULE_SIZE, reassignThreshold = 0,
                        mergeCutHeight = MERGE_CUT_HEIGHT, numericLabels = TRUE,
                        pamRespectsDendro = FALSE, deepSplit = 2, saveTOMs = FALSE, verbose = 3)

moduleColors <- labels2colors(net$colors)
module_names <- sort(unique(moduleColors))

module_sizes <- as.data.frame(table(moduleColors))
colnames(module_sizes) <- c("Module", "Gene_Count")
module_sizes <- module_sizes[order(module_sizes$Gene_Count, decreasing = TRUE), ]
write.csv(module_sizes, file.path(output_dir, "Module_Sizes.csv"), row.names = FALSE)


# -----------------------------
# 13. Panel B - gene dendrogram + modules + per-group correlation rows
#     One extra color row per treatment group: red = genes positively
#     correlated with that group, green = negatively correlated (numbers2colors).
# -----------------------------

GS_group <- cor(datExpr, trait_matrix_combined, use = "pairwise.complete.obs")
GS_group_colors <- apply(GS_group, 2, numbers2colors, signed = TRUE)
rownames(GS_group_colors) <- rownames(GS_group)

dendro_colors <- cbind(Module = moduleColors, GS_group_colors)

plot_panel_B <- function() {
  plotDendroAndColors(
    net$dendrograms[[1]], dendro_colors[net$blockGenes[[1]], ],
    groupLabels = colnames(dendro_colors),
    main = " ",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    cex.colorLabels = 1, cex.main = 1.2
  )
}
save_panel_tiff(file.path(figures_dir, "Panel_B_Dendrogram_and_Groups.tiff"),
                12, 4 + 0.25 * ncol(dendro_colors), plot_panel_B)


# -----------------------------
# 14. Module eigengenes
# -----------------------------

MEs <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)
write.csv(MEs, file.path(output_dir, "Module_Eigengenes.csv"))

ME_cor  <- cor(MEs, method = "pearson")
ME_tree <- hclust(as.dist(1 - ME_cor), method = "average")


# -----------------------------
# 15. Panel C - eigengene dendrogram + adjacency heatmap
# -----------------------------

plot_panel_C1 <- function() {
  par(mar = c(2, 4, 3, 1))
  plot(ME_tree, main = "Eigengene Dendrogram", xlab = "", sub = "", cex = 1)
}
plot_panel_C2 <- function() {
  par(mar = c(6, 6, 3, 3))
  labeledHeatmap(Matrix = ME_cor, xLabels = colnames(MEs), yLabels = colnames(MEs), ySymbols = colnames(MEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = sprintf("%.2f", ME_cor),
                 setStdMargins = FALSE, cex.text = 1, cex.lab = 1, zlim = c(-1, 1),
                 main = "Eigengene Adjacency Heatmap")
}
save_panel_tiff(file.path(figures_dir, "Panel_C_Eigengene_Network.tiff"), 11, 6,
                function() { par(mfrow = c(1, 2)); plot_panel_C1(); plot_panel_C2() })


# -----------------------------
# 16. Panel D - TOM network heatmap (random gene subsample)
# -----------------------------

set.seed(1)
tom_plot_genes <- if (ncol(datExpr) > N_GENES_TOM_PLOT) sample(colnames(datExpr), N_GENES_TOM_PLOT) else colnames(datExpr)

adjacency_sub <- adjacency(datExpr[, tom_plot_genes], power = softPower, type = NETWORK_TYPE, corFnc = CORRELATION_FUNCTION)
TOM_sub <- TOMsimilarity(adjacency_sub, TOMType = NETWORK_TYPE)
dissTOM_sub <- 1 - TOM_sub
geneTree_sub <- hclust(as.dist(dissTOM_sub), method = "average")
moduleColors_sub <- moduleColors[match(tom_plot_genes, colnames(datExpr))]

plot_matrix <- dissTOM_sub^7
diag(plot_matrix) <- NA

plot_panel_D <- function() {
  TOMplot(plot_matrix, geneTree_sub, moduleColors_sub,
          main = " ")
}
save_panel_tiff(file.path(figures_dir, "Panel_D_TOM_Heatmap.tiff"), 8, 8, plot_panel_D)

# Panels A, C, D are QC / diagnostic panels, saved standalone above.


# -----------------------------
# 17. Panel E - module-trait heatmap (factor-level traits)
# -----------------------------

moduleTraitCor <- cor(MEs, trait_matrix, use = "pairwise.complete.obs", method = "pearson")
moduleTraitPvalue <- moduleTraitCor
for (trait_name in colnames(trait_matrix)) {
  n_valid <- sum(!is.na(trait_matrix[, trait_name]))
  moduleTraitPvalue[, trait_name] <- corPvalueStudent(moduleTraitCor[, trait_name], nSamples = n_valid)
}
write.csv(moduleTraitCor, file.path(output_dir, "Module_Trait_Correlations.csv"))
write.csv(moduleTraitPvalue, file.path(output_dir, "Module_Trait_Pvalues.csv"))

textMatrixE <- sprintf("%.2f", moduleTraitCor)
dim(textMatrixE) <- dim(moduleTraitCor)

plot_panel_E <- function() {
  par(mar = c(10, 8, 4, 2))
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(trait_matrix), yLabels = colnames(MEs),
                 ySymbols = colnames(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = textMatrixE, setStdMargins = FALSE, cex.text = 0.6, cex.lab = 0.7,
                 xLabelsAngle = 90, xLabelsAdj = 1,
                 zlim = c(-1, 1), main = " ")
}
save_panel_tiff(file.path(figures_dir, "Panel_E_Module_Trait_Heatmap.tiff"), 10, 8, plot_panel_E)


# -----------------------------
# 18. Panel F - module-group heatmap, one column per treatment group
# -----------------------------

moduleGroupCor <- cor(MEs, trait_matrix_combined, use = "pairwise.complete.obs")
moduleGroupPvalue <- corPvalueStudent(moduleGroupCor, nSamples = nrow(datExpr))
write.csv(moduleGroupCor, file.path(output_dir, "Module_Group_Correlations.csv"))
write.csv(moduleGroupPvalue, file.path(output_dir, "Module_Group_Pvalues.csv"))

# Just the correlation value in each cell (no p-value line) -- keeps cells
# short so smaller, vertical column labels don't overlap the text.
textMatrixF <- sprintf("%.2f", moduleGroupCor)
dim(textMatrixF) <- dim(moduleGroupCor)

plot_panel_F <- function() {
  par(mar = c(8, 6, 1, 1.5))
  labeledHeatmap(Matrix = moduleGroupCor, xLabels = colnames(trait_matrix_combined), yLabels = colnames(MEs),
                 ySymbols = colnames(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = textMatrixF, setStdMargins = FALSE, cex.text = 0.6, cex.lab = 0.7,
                 xLabelsAngle = 90, xLabelsAdj = 1,
                 zlim = c(-1, 1), main = " ")
}
save_panel_tiff(file.path(figures_dir, "Panel_F_Module_Group_Heatmap.tiff"), 10, 8, plot_panel_F)


# -----------------------------
# 19. Gene module membership (kME) and gene significance (NA-safe)
# -----------------------------

geneModuleMembership <- cor(datExpr, MEs, use = "pairwise.complete.obs", method = "pearson")
geneModuleMembershipPvalue <- corPvalueStudent(geneModuleMembership, nSamples = nrow(datExpr))
write.csv(geneModuleMembership, file.path(output_dir, "Gene_Module_Membership_KME.csv"))
write.csv(geneModuleMembershipPvalue, file.path(output_dir, "Gene_Module_Membership_Pvalues.csv"))

geneSignificance <- matrix(NA, nrow = ncol(datExpr), ncol = ncol(trait_matrix),
                           dimnames = list(colnames(datExpr), colnames(trait_matrix)))
geneSignificancePvalue <- geneSignificance

for (trait_name in colnames(trait_matrix)) {
  trait_values <- trait_matrix[, trait_name]
  valid <- !is.na(trait_values)
  
  if (sum(valid) < 3) {
    warning("Trait '", trait_name, "' has fewer than 3 non-NA samples (", sum(valid),
            ") -- skipping, leaving NA for this trait.")
    next
  }
  
  for (gene_name in colnames(datExpr)) {
    gene_values <- datExpr[valid, gene_name]
    if (sd(gene_values) == 0 || sd(trait_values[valid]) == 0) next   # constant values -> correlation undefined
    test <- cor.test(gene_values, trait_values[valid], method = "pearson")
    geneSignificance[gene_name, trait_name] <- test$estimate
    geneSignificancePvalue[gene_name, trait_name] <- test$p.value
  }
}
write.csv(geneSignificance, file.path(output_dir, "Gene_Trait_Significance.csv"))
write.csv(geneSignificancePvalue, file.path(output_dir, "Gene_Trait_Significance_Pvalues.csv"))

# gene_info is built here (moved up from the old Table 1 section) because
# Panel I (section 23.1) needs it to pick each module's top hub genes --
# it was previously defined after Panel I, which threw
# "object 'gene_info' not found".
gene_info <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

for (ME_name in colnames(geneModuleMembership)) {
  gene_info[, ME_name] <- geneModuleMembership[, ME_name]
}

for (trait_name in colnames(geneSignificance)) {
  gene_info[, paste0("GS_", trait_name)] <- geneSignificance[, trait_name]
  gene_info[, paste0("GS_P_", trait_name)] <- geneSignificancePvalue[, trait_name]
}

write.csv(gene_info, file.path(output_dir, "Complete_Gene_WGCNA_Information.csv"), row.names = FALSE)


# -----------------------------
# 20. Focal modules (Genotype within Flower / within Leaf)
#     Reused by panels G, H, I.
# -----------------------------

best_module_for_trait <- function(trait_name) {
  cors <- moduleTraitCor[, trait_name]
  cors <- cors[!is.na(cors)]
  gsub("^ME", "", names(cors)[which.max(abs(cors))])
}

module_flower <- best_module_for_trait("Genotype_within_Flower")
module_leaf   <- best_module_for_trait("Genotype_within_Leaf")

cat("\nFocal module -- Genotype within Flower:", module_flower, "\n")
cat("Focal module -- Genotype within Leaf:  ", module_leaf, "\n")
if (module_flower == module_leaf) {
  cat("NOTE: both tissues share the same top module (", module_flower, "). Panels G/H/I for",
      "Flower and Leaf will show the same module -- that's a property of the data",
      "(this module is the strongest Genotype hit in both tissues), not a plotting bug.\n")
}


# -----------------------------
# 21. Panel G - module membership vs. gene significance
# -----------------------------

make_mm_gs_plot <- function(module, trait_name, display_label, panel_label) {
  function() {
    module_genes <- moduleColors == module
    verboseScatterplot(
      abs(geneModuleMembership[module_genes, paste0("ME", module)]),
      abs(geneSignificance[module_genes, trait_name]),
      xlab = paste("Module Membership in", module, "module"),
      ylab = paste("Gene significance for", display_label),
      main = paste0(panel_label, " ", display_label),
      cex.main = 1, cex.lab = 1, cex.axis = 1, col = module
    )
  }
}
plot_panel_G1 <- make_mm_gs_plot(module_flower, "Genotype_within_Flower", "Genotype within Flower", " ")
plot_panel_G2 <- make_mm_gs_plot(module_leaf,   "Genotype_within_Leaf",   "Genotype within Leaf",   " ")

# verboseScatterplot() can reset the device's layout internally, the same
# way plotDendroAndColors/TOMplot do for panels B/D (see the note on
# combine_panels_to_figure() in section 3). Calling both plot functions
# inside one par(mfrow = c(1, 2)) risks the second call clobbering the
# first, so each half is rendered to its own file and composited with
# magick instead -- same approach already used for panels B+F+H and the
# Panel I grid.
G1_file <- file.path(figures_dir, "Panel_G1_Flower_temp.tiff")
G2_file <- file.path(figures_dir, "Panel_G2_Leaf_temp.tiff")
save_panel_tiff(G1_file, 6, 6, plot_panel_G1)
save_panel_tiff(G2_file, 6, 6, plot_panel_G2)

combine_panels_grid(
  panel_files_matrix = list(c(G1_file, G2_file)),
  out_file = file.path(figures_dir, "Panel_G_MM_vs_GS.tiff")
)
file.remove(G1_file, G2_file)


# -----------------------------
# 22. Panel H - module expression heatmap + eigengene barplot
#
#     One module only: Flower
#     Panel size: 10" wide x 8" high (heatmap 5", barplot 3")
#
#     Widened from the original 3"x4.5": Panel H is saved standalone, at
#     the same DPI and pointsize as every other panel, so a narrower panel
#     forces the same cex to look bigger on the page than it does on the
#     wider panels (B, E, F). 10" wide keeps Panel H's text-to-panel-width
#     ratio in line with those panels, and gives the 16 treatment-group
#     labels on the barplot's x-axis enough room to not overlap.
# -----------------------------

group_breaks <- cumsum(
  rle(metadata$Combined_Treatment[sample_order])$lengths
)

group_breaks <- group_breaks[-length(group_breaks)] / nrow(metadata)

# ---------------------------------------------------------------------------
# H1. HEATMAP
# ---------------------------------------------------------------------------

plot_panel_H_heatmap <- function() {
  module <- module_flower
  
  module_genes <- colnames(datExpr)[moduleColors == module]
  
  expr_sub <- vst_matrix[
    module_genes,
    rownames(metadata),
    drop = FALSE
  ]
  
  expr_z <- t(scale(t(expr_sub)))
  
  expr_z <- expr_z[
    ,
    sample_order,
    drop = FALSE
  ]
  
  expr_z <- expr_z[
    hclust(dist(expr_z), method = "average")$order,
    ,
    drop = FALSE
  ]
  
  zmax <- max(abs(expr_z), na.rm = TRUE)
  
  par(
    mar = c(0.5, 4, 3, 1)
  )
  
  image(
    t(expr_z),
    axes = FALSE,
    col = blueWhiteRed(50),
    zlim = c(-zmax, zmax),
    main = paste0(
      module,
      " module (top ",
      nrow(expr_z),
      " genes)"
    ),
    cex.main = 1
  )
  
  abline(
    v = group_breaks,
    col = "black",
    lwd = 0.7
  )
  
  mtext(
    "Genes",
    side = 2,
    line = 1,
    cex = 0.9
  )
}

# ---------------------------------------------------------------------------
# H2. BARPLOT
# ---------------------------------------------------------------------------

plot_panel_H_barplot <- function() {
  module <- module_flower
  
  me_col <- paste0("ME", module)
  
  group_means <- tapply(
    MEs[, me_col],
    metadata$Combined_Treatment,
    mean
  )
  
  group_means <- group_means[group_order]
  
  par(
    mar = c(10, 4, 0.5, 1)
  )
  
  barplot(
    group_means,
    las = 2,
    col = module,
    ylab = "Module Eigengene",
    cex.names = 0.7,
    cex.axis = 0.8,
    cex.lab = 0.9
  )
}

# ---------------------------------------------------------------------------
# H3. SAVE PANEL H
# ---------------------------------------------------------------------------

save_panel_tiff(
  file.path(
    figures_dir,
    "Panel_H_Module_Expression_and_Eigengene.tiff"
  ),
  10,
  8,
  function() {
    layout(
      matrix(
        c(1, 2),
        nrow = 2,
        byrow = TRUE
      ),
      heights = c(5, 3)
    )
    
    plot_panel_H_heatmap()
    plot_panel_H_barplot()
  }
)


# -----------------------------
# 23. Panel I - top-5 hub gene co-expression networks
#
#     Four modules (green dropped -- too few genes for a meaningful network):
#       Panel_I_Blue.tiff
#       Panel_I_Brown.tiff
#       Panel_I_Turquoise.tiff
#       Panel_I_Yellow.tiff
#
#     Each standalone figure contains:
#       Flower network | Leaf network
#
#     All four are also combined into one 2x2 figure,
#     Panel_I_Combined_2x2.tiff (section 23.3b).
#
#     The SAME top-5 hub genes, for these same 4 modules, are used for:
#       1. Panel I network highlighting
#       2. Table 1
#
#     Hub ranking = absolute module membership (|kME|)
# -----------------------------

TOP_HUB_GENES <- 5

# Modules to show in Panel I and list in Table 1.  <- EDIT
# Module colors depend on your data, so run the script once and check Module_Sizes.csv first.
# Here green is excluded because it has too few genes for a meaningful hub network.
# Every name must be a module color that exists in your results. The 2x2 figure in section 23.3b
# is written for exactly 4 modules: edit it if you use a different number.
panel_I_modules <- c(
  "blue",
  "brown",
  "turquoise",
  "yellow"
)

# -----------------------------
# 23.1 Select top-5 hub genes
# -----------------------------

top5_hub_genes <- list()

for (module in panel_I_modules) {
  ME_name <- paste0("ME", module)
  
  module_genes <- gene_info[
    gene_info$Module == module,
    ,
    drop = FALSE
  ]
  
  if (!ME_name %in% colnames(module_genes)) {
    warning("Module eigengene not found for module: ", module)
    next
  }
  
  # Absolute module membership = kME
  module_genes$KME <- abs(module_genes[[ME_name]])
  
  # Highest kME first
  module_genes <- module_genes[
    order(module_genes$KME, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  # EXACTLY the top 5 genes
  module_genes <- head(module_genes, TOP_HUB_GENES)
  
  # Rank 1-5
  module_genes$Rank <- seq_len(nrow(module_genes))
  
  top5_hub_genes[[module]] <- module_genes
}

# Combine all modules into one table
top5_hub_table <- do.call(rbind, top5_hub_genes)
rownames(top5_hub_table) <- NULL

# Save the exact gene list used in Panel I and Table 1
write.csv(
  top5_hub_table,
  file.path(output_dir, "Top_5_Hub_Genes_Per_Module.csv"),
  row.names = FALSE
)

# -----------------------------
# 23.2 Network function
# -----------------------------

make_hub_network <- function(module, tissue_code, tissue_label) {
  # Samples from the selected tissue
  sample_ids <- rownames(metadata)[metadata$Tissue == tissue_code]
  
  if (length(sample_ids) < 3) {
    plot.new()
    title(paste(module, tissue_label, "- insufficient samples"))
    return()
  }
  
  # Genes belonging to this WGCNA module
  module_genes <- colnames(datExpr)[moduleColors == module]
  
  if (length(module_genes) < 2) {
    plot.new()
    title(paste(module, tissue_label, "- insufficient genes"))
    return()
  }
  
  # Tissue-specific expression
  expr_module <- datExpr[
    sample_ids,
    module_genes,
    drop = FALSE
  ]
  
  # ------------------------------------------------------------
  # EXACT SAME TOP-5 HUB GENES USED IN TABLE 1
  # ------------------------------------------------------------
  hub_genes <- top5_hub_genes[[module]]$Gene
  
  # Keep only genes available in this tissue/network
  hub_genes <- hub_genes[
    hub_genes %in% colnames(expr_module)
  ]
  
  # ------------------------------------------------------------
  # Calculate adjacency and TOM
  # ------------------------------------------------------------
  adj <- adjacency(
    expr_module,
    power = softPower,
    type = NETWORK_TYPE,
    corFnc = CORRELATION_FUNCTION
  )
  
  tom <- TOMsimilarity(
    adj,
    TOMType = NETWORK_TYPE
  )
  
  dimnames(tom) <- list(
    colnames(expr_module),
    colnames(expr_module)
  )
  
  # Remove weak edges
  tom[tom < EDGE_WEIGHT_THRESHOLD] <- 0
  diag(tom) <- 0
  
  # ------------------------------------------------------------
  # Keep hub genes + genes connected to hubs
  # ------------------------------------------------------------
  connected_genes <- unique(
    c(
      hub_genes,
      unlist(
        lapply(
          hub_genes,
          function(g) {
            names(tom[g, ])[tom[g, ] > 0]
          }
        )
      )
    )
  )
  
  connected_genes <- connected_genes[
    connected_genes %in% colnames(expr_module)
  ]
  
  if (length(connected_genes) < 2) {
    plot.new()
    title(paste(module, tissue_label, "- no network"))
    return()
  }
  
  tom_sub <- tom[
    connected_genes,
    connected_genes,
    drop = FALSE
  ]
  
  # ------------------------------------------------------------
  # Build graph
  # ------------------------------------------------------------
  g <- graph_from_adjacency_matrix(
    tom_sub,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  # Remove isolated nodes
  g <- delete_vertices(
    g,
    degree(g) == 0
  )
  
  # ------------------------------------------------------------
  # Node appearance
  # ------------------------------------------------------------
  V(g)$color <- module
  V(g)$size <- 4
  V(g)$label <- NA
  
  # Highlight EXACT top-5 hub genes
  hub_nodes <- V(g)$name %in% hub_genes
  
  V(g)$color[hub_nodes] <- "red"
  V(g)$size[hub_nodes] <- 7
  V(g)$label[hub_nodes] <- V(g)$name[hub_nodes]
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  par(
    mar = c(1, 1, 3, 1)
  )
  
  set.seed(123)
  
  plot(
    g,
    layout = layout_with_fr(g),
    vertex.label.color = "black",
    vertex.label.cex = 0.65,
    vertex.frame.color = V(g)$color,
    edge.width = E(g)$weight * 3,
    main = paste0(
      tissue_label,
      " - ",
      module,
      " module"
    )
  )
}

# -----------------------------
# 23.3 Create one figure for each module
#
#     Flower = F
#     Leaf   = L
#
#     Each figure = 6.27 x 3 inches
# -----------------------------

for (module in panel_I_modules) {
  plot_flower <- function() {
    make_hub_network(
      module = module,
      tissue_code = "F",
      tissue_label = "Flower"
    )
  }
  
  plot_leaf <- function() {
    make_hub_network(
      module = module,
      tissue_code = "L",
      tissue_label = "Leaf"
    )
  }
  
  output_file <- file.path(
    figures_dir,
    paste0(
      "Panel_I_",
      tools::toTitleCase(module),
      ".tiff"
    )
  )
  
  save_panel_tiff(
    output_file,
    6.27,
    3,
    function() {
      par(
        mfrow = c(1, 2)
      )
      
      plot_flower()
      plot_leaf()
    }
  )
}

# -----------------------------
# 23.3b Combine the 4 module figures into one 2x2 figure
# -----------------------------

panel_I_files <- setNames(
  file.path(
    figures_dir,
    paste0("Panel_I_", tools::toTitleCase(panel_I_modules), ".tiff")
  ),
  panel_I_modules
)

combine_panels_grid(
  panel_files_matrix = list(
    c(panel_I_files[["blue"]],      panel_I_files[["brown"]]),
    c(panel_I_files[["turquoise"]], panel_I_files[["yellow"]])
  ),
  out_file = file.path(figures_dir, "Panel_I_Combined_2x2.tiff")
)

# -----------------------------
# 23.4 Optional: list all panel-I figures created
# -----------------------------

cat("\nPanel I figures created:\n")

for (module in panel_I_modules) {
  cat(
    "  ",
    file.path(
      figures_dir,
      paste0(
        "Panel_I_",
        tools::toTitleCase(module),
        ".tiff"
      )
    ),
    "\n"
  )
}

cat("  ", file.path(figures_dir, "Panel_I_Combined_2x2.tiff"), "(combined)\n")


# -----------------------------
# 24. Table 1 - top-5 hub genes with gene annotation
#
#     IMPORTANT:
#     Table 1 uses EXACTLY the same top-5 genes selected above for Panel I.
#
#     Therefore:
#       Panel I hub genes = Table 1 hub genes
# -----------------------------

# -----------------------------
# 24.2 Save all kME-threshold hub genes as supplementary data
#
#     This section is independent of the top-5 selection.
#     It preserves all genes satisfying KME_THRESHOLD.
# -----------------------------

hub_dir <- file.path(
  output_dir,
  "Hub_Genes"
)

if (!dir.exists(hub_dir)) {
  dir.create(
    hub_dir,
    recursive = TRUE
  )
}

hub_gene_tables <- list()

for (module in module_names) {
  ME_name <- paste0("ME", module)
  
  module_genes <- gene_info[
    gene_info$Module == module,
    ,
    drop = FALSE
  ]
  
  if (ME_name %in% colnames(module_genes)) {
    module_genes$KME <-
      abs(module_genes[[ME_name]])
    
    module_genes <- module_genes[
      order(
        module_genes$KME,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ]
    
    hub_genes <- module_genes[
      module_genes$KME >= KME_THRESHOLD,
      ,
      drop = FALSE
    ]
    
    hub_gene_tables[[module]] <- hub_genes
    
    write.csv(
      hub_genes,
      file.path(
        hub_dir,
        paste0(
          "Hub_Genes_",
          module,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
}

# -----------------------------
# 24.3 Table 1 = exact top-5 genes used in Panel I
# -----------------------------

Table1_hub_genes <- top5_hub_table

# -----------------------------
# 24.4 GFF annotation
# -----------------------------

if (nrow(Table1_hub_genes) == 0) {
  warning(
    "No top hub genes were selected. ",
    "Table 1 was not created."
  )
} else if (file.exists(gff_file)) {
  get_field <- function(x, key) {
    pattern <- paste0(
      key,
      "=([^;]+)"
    )
    
    r <- regexpr(
      pattern,
      x
    )
    
    m <- regmatches(
      x,
      r
    )
    
    result <- rep(
      NA_character_,
      length(x)
    )
    
    result[r > 0] <- m
    
    sub(
      paste0("^", key, "="),
      "",
      result
    )
  }
  
  # Read mRNA records
  mrna_lines <- grep(
    "\tmRNA\t",
    readLines(gff_file),
    value = TRUE
  )
  
  if (length(mrna_lines) == 0) {
    warning(
      "No '\\tmRNA\\t' lines found in GFF file: ",
      gff_file,
      "\nTable 1 annotation was not created."
    )
  } else {
    attrs <- vapply(
      strsplit(
        mrna_lines,
        "\t"
      ),
      `[`,
      character(1),
      9
    )
    
    # Gene ID
    gene_id <- sub(
      "^gene:",
      "",
      get_field(
        attrs,
        "Parent"
      )
    )
    
    # Functional description
    description <- get_field(
      attrs,
      "Note"
    )
    
    annotation <- data.frame(
      Gene_ID = gene_id,
      Description = description,
      stringsAsFactors = FALSE
    )
    
    annotation <- annotation[
      !is.na(annotation$Gene_ID) &
        !duplicated(annotation$Gene_ID),
      ,
      drop = FALSE
    ]
    
    write.csv(
      annotation,
      file.path(
        output_dir,
        "ITAG4.0_Gene_Annotation.csv"
      ),
      row.names = FALSE
    )
    
    ###########################################################################
    # MERGE TOP-5 HUB GENES WITH ANNOTATION
    ###########################################################################
    
    Table1 <- merge(
      Table1_hub_genes[
        ,
        c(
          "Module",
          "Gene",
          "KME",
          "Rank"
        ),
        drop = FALSE
      ],
      annotation,
      by.x = "Gene",
      by.y = "Gene_ID",
      all.x = TRUE
    )
    
    # Sort by module and rank
    Table1 <- Table1[
      order(
        Table1$Module,
        Table1$Rank
      ),
      ,
      drop = FALSE
    ]
    
    ###########################################################################
    # SAVE TABLE 1
    ###########################################################################
    
    write.csv(
      Table1,
      file.path(
        output_dir,
        "Table_1_Top_5_Hub_Genes.csv"
      ),
      row.names = FALSE
    )
    
    ###########################################################################
    # ALSO KEEP THE ORIGINAL TABLE NAME
    ###########################################################################
    
    write.csv(
      Table1,
      file.path(
        output_dir,
        "Hub_Gene_Annotation_Table.csv"
      ),
      row.names = FALSE
    )
    
    ###########################################################################
    # ANNOTATION MATCH RATE
    ###########################################################################
    
    matched_pct <- round(
      100 *
        mean(
          !is.na(Table1$Description)
        ),
      1
    )
    
    cat(
      "\nTable 1 annotation match rate: ",
      matched_pct,
      "%\n",
      sep = ""
    )
  }
} else {
  warning(
    "GFF file not found at: ",
    gff_file,
    "\nTable 1 annotation was not created."
  )
}

# -----------------------------
# 24.5 Print the exact genes used in Panel I and Table 1
# -----------------------------

cat("\n============================================================\n")
cat("TOP-5 HUB GENES USED IN BOTH PANEL I AND TABLE 1\n")
cat("============================================================\n")

for (module in panel_I_modules) {
  cat(
    "\n",
    toupper(module),
    " module:\n",
    sep = ""
  )
  
  if (!is.null(top5_hub_genes[[module]])) {
    print(
      top5_hub_genes[[module]][
        ,
        c(
          "Rank",
          "Gene",
          "KME"
        ),
        drop = FALSE
      ]
    )
  } else {
    cat("No genes available.\n")
  }
}

cat("\n============================================================\n")
cat("Panel I and Table 1 use the EXACT SAME top-5 hub genes.\n")
cat("============================================================\n")


# -----------------------------
# 25. Module size bar plot (supplementary)
# -----------------------------

module_sizes$Module <- factor(module_sizes$Module, levels = module_sizes$Module)
p_module_size <- ggplot(module_sizes, aes(x = Module, y = Gene_Count)) +
  geom_col() +
  theme_classic(base_size = 14) +
  labs(title = "WGCNA Module Sizes", x = "Module", y = "Number of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11), axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))
ggsave(file.path(output_dir, "07_Module_Size.pdf"), p_module_size, width = 9, height = 6)


# -----------------------------
# 26. Save objects
# -----------------------------

saveRDS(net, file.path(output_dir, "WGCNA_Network.rds"))
saveRDS(datExpr, file.path(output_dir, "WGCNA_Expression_Matrix.rds"))
saveRDS(MEs, file.path(output_dir, "Module_Eigengenes.rds"))

save(counts, counts_filtered, vst_matrix, datExpr, metadata,
     trait_matrix, trait_matrix_combined, net, moduleColors, MEs,
     moduleTraitCor, moduleTraitPvalue, moduleGroupCor, moduleGroupPvalue,
     geneModuleMembership, geneModuleMembershipPvalue,
     geneSignificance, geneSignificancePvalue, gene_info, module_sizes, softPower,
     file = file.path(output_dir, "Complete_WGCNA_Analysis.RData"))


# -----------------------------
# 27. Final summary
# -----------------------------

cat("\nWGCNA finished.\n")
cat("Genes used:", ncol(datExpr), " Samples:", nrow(datExpr), " Modules:", length(module_names), " Power:", softPower, "\n")
cat("Focal modules -- Genotype within Flower:", module_flower, " | Genotype within Leaf:", module_leaf, "\n")
cat("Figures:\n ", normalizePath(figures_dir), "\n")
cat("Panels A-H saved standalone as Panel_*.tiff (no multi-panel composite)\n")
cat("Panel I: 4 standalone module figures (green dropped) + Panel_I_Combined_2x2.tiff\n")
cat("Table 1 equivalent: Hub_Gene_Annotation_Table.csv, hub genes match panel I exactly (requires gff_file to exist)\n")