##Bulk RNA Seq Data Analysis Workflow by Md Jahid Hasan Jone##

library(rtracklayer)
library(dplyr)

# Import GFF from reference folder
gff <- import("/.../.../Bulk_RNA_seq/2_References/gene_annotation.gff")

# Keep mRNA features only
mrna <- gff[gff$type == "mRNA"]


# Extract transcript and gene IDs
# Parent is a list column in rtracklayer, so convert it to a plain character vector

tx2gene <- data.frame(
  TXNAME = as.character(mrna$ID),
  GENEID = as.character(unlist(mrna$Parent))
)

# Remove "mRNA:" and "gene:" prefixes if present (e.g. "mRNA:Solyc01g005000.3.1"),
# so IDs match the fasta headers / quant.sf names. Does nothing if there is no prefix.
tx2gene$TXNAME <- sub("^mRNA:", "", tx2gene$TXNAME)
tx2gene$GENEID <- sub("^gene:", "", tx2gene$GENEID)

# Remove duplicates
tx2gene <- unique(tx2gene)

head(tx2gene)

write.csv(
  tx2gene,
  "/.../.../Bulk_RNA_seq/2_References/tx2gene.csv",
  row.names = FALSE
)



#####Validate#####

library(readr)

tx2gene <- read_csv("/.../.../Bulk_RNA_seq/2_References/tx2gene.csv")
quant <- read.delim(
  "/.../.../Bulk_RNA_seq/5_Results/2_Salmon/.../quant.sf"    #quant.sf file location for any sample
)
sum(quant$Name %in% tx2gene$TXNAME)

# Fraction of transcripts matched (should be close to 1)
mean(quant$Name %in% tx2gene$TXNAME)