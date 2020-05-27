#!/bin/bash -l
#SBATCH
#SBATCH --job-name=compute_pseudobulk
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=1:0:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4

module load gcc/5.5.0
module load R/3.6.1
Rscript ../src/processing_pseudobulk.R $1 > ../log/processing_pseudobulk/$1.log 
