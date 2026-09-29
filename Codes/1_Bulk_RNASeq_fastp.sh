#!/bin/bash
#BSUB -n 4
#BSUB -W 24:00
#BSUB -J fastp
#BSUB -o stdout.%J
#BSUB -e stderr.%J

INPUTDIR=/.../.../Bulk_RNA_seq/3_Raw_Reads
OUTPUTDIR=/.../.../Bulk_RNA_seq/4_Trimmed_Reads
RESULTDIR=/.../.../Bulk_RNA_seq/5_Results/1_fastp

mkdir -p $INPUTDIR
mkdir -p $OUTPUTDIR
mkdir -p $RESULTDIR

module load conda
source ~/.bashrc
conda activate /.../.../usrapps/group/gatk_rnaseq   # conda location


for file1 in "$INPUTDIR"/*_1.fq.gz
do
    # Extract sample name
    # Example:
    #   2F_T0_48h_R1_1.fq.gz
    #   ↓
    #   sample = 2F_T0_48h_R1
    sample=$(basename "$file1" "_1.fq.gz")

    filename1=${sample}_1.fq.gz
    filename2=${sample}_2.fq.gz

    if [ -f $INPUTDIR/$filename1 ] && [ -f $INPUTDIR/$filename2 ];

	then
	# -I, -O	: paired end only (remove both lines for single end)
        fastp \
        -i $INPUTDIR/$filename1 \
        -I $INPUTDIR/$filename2 \
        -o $OUTPUTDIR/TMD_$filename1 \
        -O $OUTPUTDIR/TMD_$filename2 \
        -h $RESULTDIR/${sample}.html \
        -j $RESULTDIR/${sample}.json
    fi
done

conda deactivate
