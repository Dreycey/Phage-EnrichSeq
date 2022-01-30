#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=sknl
#SBATCH --output=sample-%j.out

module purge

# update log file
LOG='Benchmarking.log'
echo "== This is the scripting step! ==" >> ${LOG}
echo "Benchmaarking " >> ${LOG}

# load modules
echo "initiating conda environment and dependencies" >> ${LOG};
source ~/.bash_profile;
conda activate enrichseq;

# sttart benchmarking
echo "Running Testing now" >> ${LOG}
bash run_benchmarking.sh;
# send note that it finished
echo "Benchmarking all tools has finished"
