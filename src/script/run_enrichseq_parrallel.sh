#!/bin/bash

#SBATCH --time 12:00:00
#SBATCH --partition smem
#SBATCH --ntasks=10
#SBATCH --job-name EnrichSeq_T3_k28_parrallel
#SBATCH --output enrichbench_T3_k28_parrallel.out

module purge
module load anaconda 
conda activate enrichseq
module load gnu_parallel


#######
# Choose test type
#######
#bash RUNTESTS_3_onebyone.txt;
parallel < RUNTESTS_3.txt;
