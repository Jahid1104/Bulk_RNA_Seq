library(rtracklayer)
library(dplyr)

# Import GFF
#gff <- import("annotation.gff") from reference folder
gff <- import(file.choose())

# Keep mRNA features only
mrna <- gff[gff$type == "mRNA"]


# Extract transcript and gene IDs

tx2gene <- data.frame(
  TXNAME = mrna$ID,
  GENEID = mrna$Parent
)

# Remove duplicates
tx2gene <- unique(tx2gene)

head(tx2gene)

write.csv(
  tx2gene,
  "tx2gene.csv",
  row.names = FALSE
)



#####Validate#####

library(readr)

tx2gene <- read_csv("tx2gene.csv")
quant <- read.delim(
  "1F_T1_48h_R2/quant.sf"
)
sum(quant$Name %in% tx2gene$TXNAME)
