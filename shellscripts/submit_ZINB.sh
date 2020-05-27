#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_ZINB
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:10:0
#SBATCH --partition=express
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

for ngenes in 10
do
	for learning_rate in 0.0001 0.00001
	do 
		for max_epochs in 80000 100000 150000  
			do
        			sbatch compute_ZINB.sh $ngenes $learning_rate $max_epochs
			done
	done
done
