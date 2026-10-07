#!/bin/bash
#BSUB -W 20
#BSUB -n 1
##BSUB -R "rusage[mem=??GB]"   #Specify maximum memory required
##BSUB -x                      #Use exclusive only if necessary
#BSUB -J WGCNA
#BSUB -o out.%J
#BSUB -e err.%J

module load R
Rscript /.../.../Bulk_RNA_seq/1_Codes/8_WGCNA.R
