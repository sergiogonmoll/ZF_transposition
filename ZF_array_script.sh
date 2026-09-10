#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=16-89
#SBATCH --job-name=TEF
#SBATCH --mem-per-cpu=250GB
#SBATCH --partition=regular


# Make all the directories needed
mkdir ${SLURM_ARRAY_TASK_ID}RetroSeq

module load BEDTools
module load SAMtools/0.1.19-GCC-10.3.0 
module load Exonerate

# Run RetroSeq discovery #
perl retroseq.pl -discover -bam Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_all.bam -output ${SLURM_ARRAY_TASK_ID}RetroSeq/SRR304140${SLURM_ARRAY_TASK_ID}.tab -refTEs refTEs.tab -eref TEfasta.tab -align

# Run RetroSeq call #
perl retroseq.pl -call -bam Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_all.bam -input ${SLURM_ARRAY_TASK_ID}RetroSeq/SRR304140${SLURM_ARRAY_TASK_ID}.tab  -ref GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -output ${SLURM_ARRAY_TASK_ID}RetroSeq/SRR304140${SLURM_ARRAY_TASK_ID}_TE.vcf  -reads 10 -depth 400 -soft
#gzip ${SLURM_ARRAY_TASK_ID}RetroSeq/SRR304140${SLURM_ARRAY_TASK_ID}_TE.vcf


#gunzip focal.vcf.gz
#gunzip dam.vcf.gz
#gunzip sire.vcf.gz
#module load R
#Rscript VCF_RetroSeq.R
