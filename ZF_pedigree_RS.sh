#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-45
#SBATCH --job-name=TEF
#SBATCH --mem-per-cpu=10GB
#SBATCH --partition=regular

array_number=${1} # Use this as a way to read the correct file

dam=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${array_number}| awk '{print $1}')

sire=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${array_number}| awk '{print $2}')

focal=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${array_number}| awk '{print $3}')


# Make all the directories needed
mkdir TEFLON/Results/

module load R

Rscript RetroSeq_output.R ${dam} ${sire} ${focal}

#cp ${focal}private_TEs.csv Results/
#cp ${focal}private_TEs_detail.csv Results/
cp ${focal}TEFLoN_TEs_detail.csv TEFLON/Results/
