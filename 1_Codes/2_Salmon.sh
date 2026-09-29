#!/bin/bash
#BSUB -n 4
#BSUB -W 48:00
#BSUB -R select[avx2]
##BSUB -R "rusage[mem=64GB]"
#BSUB -J salmon
#BSUB -o out.%J
#BSUB -e err.%J

seq_path=/.../.../Bulk_RNA_seq/2_References
transcriptome=$seq_path/ITAG4.0_cDNA.fasta      # update with your reference transcriptome fasta file
salmon_tmt_index=$seq_path/salmon_tmt_index

INPUTDIR=/.../.../Bulk_RNA_seq/4_Trimmed_Reads
RESULTDIR=/.../.../Bulk_RNA_seq/5_Results/2_Salmon

mkdir -p $salmon_tmt_index
mkdir -p $RESULTDIR

module load conda
source ~/.bashrc
conda activate /.../.../usrapps/group/gatk_rnaseq   # conda location


# Create transcriptome index (only if not already built).
if [ ! -f $salmon_tmt_index/versionInfo.json ]; then
    salmon index -t $transcriptome -i $salmon_tmt_index -p 4
fi


for file1 in "$INPUTDIR"/TMD_*_1.fq.gz
do
    sample=$(basename "$file1" "_1.fq.gz")
    sample=${sample#TMD_}

    filename1=TMD_${sample}_1.fq.gz
    filename2=TMD_${sample}_2.fq.gz

    if [ -f $INPUTDIR/$filename1 ] && [ -f $INPUTDIR/$filename2 ];
    then
        echo "Quantifying: $sample"
		# -1 and -2 : paired end. For single end, use -r $INPUTDIR/$filename1 instead.
        salmon quant \
        -i $salmon_tmt_index \
        --libType A \
        -1 $INPUTDIR/$filename1 \
        -2 $INPUTDIR/$filename2 \
        -p 4 \
        --seqBias \
        -o $RESULTDIR/$sample
    fi
done

conda deactivate

