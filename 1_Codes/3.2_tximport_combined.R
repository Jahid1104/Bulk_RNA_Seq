# tximport_combined.R
# tximport + expression plots for ALL Salmon outputs combined

# -----------------------------
# Packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

if (!requireNamespace("tximport", quietly = TRUE)) 
  BiocManager::install("tximport", update = FALSE, ask = FALSE)

pkgs <- c("readr", "dplyr", "tidyr", "tibble", "ggplot2", "stringr", "DESeq2")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) 
  install.packages(to_install)


library(tximport)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)
library(DESeq2)



# -----------------------------
# User paths
# -----------------------------
setwd("/.../.../Bulk_RNA_seq")

salmon_dir <- "5_Results/2_Salmon"
tx2gene_file  <- "2_References/tx2gene.csv"
output_dir    <- "5_Results/3_tximport"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper function: extract tissue (F/L) from a sample name
# e.g. "2F_T0_48h_R1" -> "F", "2L_T1_72h_R3" -> "L"
# -----------------------------
extract_tissue <- function(sample_names) {
  str_match(sample_names, "^\\d+([FL])")[, 2]
}

# -----------------------------
# Helper function: natural sort of sample names
# Pads every run of digits so "10F" sorts after "2F" (not before it),
# and works for any sample name shape - not just single letters/numbers.
# -----------------------------
sort_samples <- function(x) {
  keys <- str_replace_all(x, "\\d+", function(m) formatC(as.integer(m), width = 10, flag = "0"))
  x[order(keys)]
}

# -----------------------------
# Helper function: get quant.sf files from one or more Salmon result roots
# -----------------------------
get_quant_files <- function(root_dirs) {
  files <- unlist(lapply(root_dirs, function(d) {
    list.files(d, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
  }))
  sample_names <- basename(dirname(files))
  names(files) <- sample_names
  
  # Guard against the same sample name showing up in both the PE and SE
  # folders (e.g. a sample that got reprocessed). tximport requires unique
  # names, and silently keeping only one copy could quietly drop data, so
  # fail loudly instead and say exactly which names collided.
  dupes <- unique(sample_names[duplicated(sample_names)])
  if (length(dupes) > 0) {
    stop(
      "Duplicate sample name(s) found across Salmon_PE and Salmon_SE: ",
      paste(dupes, collapse = ", "),
      ". Each sample should have a quant.sf in only one of the two folders."
    )
  }
  
  sample_names <- sort_samples(sample_names)
  files <- files[sample_names]
  files
}

# -----------------------------
# Helper function: make plots and save TIFFs
# -----------------------------
make_tximport_plots <- function(txi, run_label, outdir) {
  sample_order <- colnames(txi$counts)
  tissue_lookup <- setNames(extract_tissue(sample_order), sample_order)
  
  counts_long <- as.data.frame(txi$counts) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
    mutate(sample = factor(sample, levels = sample_order),
           tissue = tissue_lookup[as.character(sample)],
           log_count = log2(count + 1))
  
  tpm_long <- as.data.frame(txi$abundance) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>%
    mutate(sample = factor(sample, levels = sample_order),
           tissue = tissue_lookup[as.character(sample)],
           log_tpm = log2(tpm + 1))
  
  sample_summary <- data.frame(
    sample = factor(sample_order, levels = sample_order),
    tissue = tissue_lookup[sample_order],
    total_counts = colSums(txi$counts),
    total_tpm = colSums(txi$abundance)
  )
  
  p1 <- ggplot(counts_long, aes(x = sample, y = log_count)) +
    geom_boxplot(fill = "steelblue", outlier.size = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "log2(gene counts + 1)", title = " ")
  
  # Colored by tissue (F/L) rather than by individual sample - with 70+
  # samples now combined, per-sample coloring is too crowded to read.
  p2 <- ggplot(counts_long, aes(x = log_count, color = tissue, fill = tissue)) +
    geom_density(alpha = 0.12) +
    theme_bw() +
    labs(x = "log2(gene counts + 1)", y = "Density", title = " ")
  
  p3 <- ggplot(tpm_long, aes(x = sample, y = log_tpm)) +
    geom_boxplot(fill = "darkgreen", outlier.size = 0.2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "log2(TPM + 1)", title = " ")
  
  p4 <- ggplot(tpm_long, aes(x = log_tpm, color = tissue, fill = tissue)) +
    geom_density(alpha = 0.12) +
    theme_bw() +
    labs(x = "log2(TPM + 1)", y = "Density", title = " ")
  
  p5 <- ggplot(sample_summary, aes(x = sample, y = total_counts, fill = tissue)) +
    geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Sample", y = "Total gene counts", title = " ")
  
  # PCA
  coldata <- data.frame(row.names = sample_order, tissue = tissue_lookup[sample_order])
  dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ 1)
  
  # With 72 samples spanning genotype, tissue, temperature, and time, the
  # gene-wise dispersion trend can fail to fit a parametric curve, which
  # makes vst() error out and stops the function before p6 (and the six
  # save_tiff() calls below it) ever run. Fall back to a local fit instead
  # of letting the whole script die partway through.
  vsd <- tryCatch(
    vst(dds, blind = TRUE),
    error = function(e) {
      message(
        run_label, ": parametric vst() fit failed (",
        conditionMessage(e), "), falling back to fitType = 'local'."
      )
      varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
    }
  )
  pca <- prcomp(t(assay(vsd)))
  
  pca_df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) %>%
    mutate(sample = factor(sample, levels = sample_order),
           tissue = tissue_lookup[as.character(sample)])
  
  p6 <- ggplot(pca_df, aes(PC1, PC2, label = sample, color = tissue)) +
    geom_point(size = 3) +
    geom_text(nudge_y = 0.02, size = 3, show.legend = FALSE) +
    theme_bw() +
    labs(title = " ")
  
  save_tiff <- function(plot_obj, filename, width = 10, height = 6) {
    ggsave(
      filename = file.path(outdir, filename),
      plot = plot_obj,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      compression = "lzw"
    )
  }
  
  save_tiff(p1, paste0(run_label, "_gene_count_distribution.tiff"))
  save_tiff(p2, paste0(run_label, "_gene_count_density.tiff"))
  save_tiff(p3, paste0(run_label, "_tpm_distribution.tiff"))
  save_tiff(p4, paste0(run_label, "_tpm_density.tiff"))
  save_tiff(p5, paste0(run_label, "_total_counts_per_sample.tiff"))
  save_tiff(p6, paste0(run_label, "_pca_samples.tiff"), width = 8, height = 6)
  
  invisible(list(
    counts_long = counts_long,
    tpm_long = tpm_long,
    sample_summary = sample_summary
  ))
}

# -----------------------------
# Read tx2gene
# -----------------------------
tx2gene <- read_csv(tx2gene_file, col_names = c("TXNAME", "GENEID"), show_col_types = FALSE)

# -----------------------------
# Get files (PE + SE combined - all 72 samples together)
# -----------------------------
combined_files <- get_quant_files(c(salmon_pe_dir, salmon_se_dir))

# -----------------------------
# tximport
# -----------------------------
txi_combined <- tximport(
  combined_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

saveRDS(
  txi_combined,
  file.path(output_dir, "txi.rds")
)

# -----------------------------
# Write output tables
# -----------------------------
write.csv(txi_combined$counts, file.path(output_dir, "gene_counts.csv"))
write.csv(txi_combined$abundance, file.path(output_dir, "TPM.csv"))

# -----------------------------
# Create plots
# -----------------------------
combined_results <- make_tximport_plots(txi_combined, "combined", output_dir)

# -----------------------------
# Save sample summary
# -----------------------------
write.csv(combined_results$sample_summary, file.path(output_dir, "sample_summary.csv"), row.names = FALSE)

message("Done. Outputs saved to: ", output_dir)
