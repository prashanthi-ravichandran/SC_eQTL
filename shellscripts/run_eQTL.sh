#!/bin/bash -l
#SBATCH
#SBATCH --job-name=eQTL
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:15:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=3

module load gcc/5.5.0
module load R/3.6.1
Rscript ../src/call_eQTLs.R $1 $2 > ../log/calling_eQTLs_simple/$1/chr$2.log
