#!/bin/bash
#SBATCH --job-name=BWA_align
#SBATCH --partition=regular
#SBATCH --cpus-per-task=6
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --output=BWA_SW-%j.log
#SBATCH --array=17-89
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.a.gonzalez.mollinedo@rug.nl

module load BWA
module load SAMtools

# make directories
# only need to run once
#mkdir Aligned_ZF/
#mkdir /scratch/p309374/SW_genomics/Aligned/${data_name}/

# index file to be used by bwa
# only need to run once
#bwa index /home3/p309374/Zebra_Finch/Ref_genome/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz

# align paired reads using bwa mem and output as bam file using samtools
bwa mem -t 6 \
        /home3/p309374/Zebra_Finch/Ref_genome/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz \
        Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_paired_R1.fastq.gz \
        Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_paired_R2.fastq.gz | \
        samtools sort -o Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_paired.bam

# combine single end reads into one file
# align unpaired reads using bwa mem and output as bam file using samtools
# merge paired and unpaired alignments



zcat  Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_unpaired_R1.fastq.gz Trimmed_ZF/ Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_unpaired_R2.fastq.gz > Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_unpaired_both.fastq.gz
bwa mem -t 6 \
		/home3/p309374/Zebra_Finch/Ref_genome/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz \
		Trimmed_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_trimmed_unpaired_both.fastq.gz| \
		samtools sort -o Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_unpaired.bam
		samtools merge -@ 6 Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_all.bam Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_paired.bam Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_unpaired.bam

samtools index -@ 6 Aligned_ZF/SRR304140${SLURM_ARRAY_TASK_ID}_all.bam
