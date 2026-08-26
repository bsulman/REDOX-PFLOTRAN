#!/bin/bash

#SBATCH --account=ccsi
#SBATCH -p batch 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 6:00:00
#SBATCH -o Arctic_Fe_output/logfiles/Arctic_out_%A-%a.txt
#SBATCH -e Arctic_Fe_output/logfiles/Arctic_error_%A-%a.txt

# Run an array of jobs with numbers 0-4 that will be set in SLURM_ARRAY_TASK_ID
# To limit number of simultaneous jobs use %N, e.g. 0-32%8 for max 8 jobs to run at once
#SBATCH --array=1-8
## SBATCH --ntasks-per-node=2

# Nodes appear to have 32 or 36 CPUs and it looks like it will schedule up to that many on a single node from the array
# 2 sockets/16 or 18 cores/1 thread per node (%z)
# sinfo -p batch -o "%.6D %.11T %.4c %.8z %.6m %.8d %.6w %.8f %20E"

source ~/.bashrc

module load anaconda3
conda activate py39f

python Arctic_Fe_CH4.py -n ${SLURM_ARRAY_TASK_ID} -N ${SLURM_ARRAY_TASK_MAX}
# python Arctic_Fe_CH4_paramtest.py -n ${SLURM_ARRAY_TASK_ID} -N ${SLURM_ARRAY_TASK_MAX}