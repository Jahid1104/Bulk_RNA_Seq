##Bulk RNA Seq Data Analysis Workflow by Md Jahid Hasan Jone##

library(tximport)
library(readr)
library(stringr)

# Paths
salmon_dir   <- "/.../.../Bulk_RNA_seq/5_Results/2_Salmon"
tx2gene_file <- "/.../.../Bulk_RNA_seq/2_References/tx2gene.csv"
output_dir   <- "/.../.../Bulk_RNA_seq/5_Results/3_tximport"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# Natural sort of sample names
# Pads every run of digits so "10F" sorts after "2F" (not before it)
sort_samples <- function(x) {
  keys <- str_replace_all(x, "\\d+", function(m) formatC(as.integer(m), width = 10, flag = "0"))
  x[order(keys)]
}


# Get quant.sf files (one per sample folder) and name them by sample
files <- list.files(salmon_dir, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
files <- files[sort_samples(names(files))]


# Read tx2gene (first two columns: TXNAME, GENEID)
tx2gene <- read_csv(tx2gene_file, show_col_types = FALSE)


# tximport
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)


# Write outputs
saveRDS(txi, file.path(output_dir, "txi.rds"))
write.csv(txi$counts, file.path(output_dir, "gene_counts.csv"))
write.csv(txi$abundance, file.path(output_dir, "TPM.csv"))

message("Done. ", length(files), " samples. Outputs saved to: ", output_dir)