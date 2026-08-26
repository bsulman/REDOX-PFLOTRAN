#!/bin/bash

#SBATCH --account=ccsi
#SBATCH -p batch 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH -t 8:00:00
#SBATCH -o /lustre/or-scratch/cades-ccsi/b0u/Mn_output/logfiles/Mn_out_%A-%a.txt
#SBATCH -e /lustre/or-scratch/cades-ccsi/b0u/Mn_output/logfiles/Mn_error_%A-%a.txt

# Run an array of jobs with numbers 0-4 that will be set in ${SLURM_ARRAY_TASK_ID}
# To limit number of simultaneous jobs use %N, e.g. 0-32%8 for max 8 jobs to run at once
#SBATCH --array=1-32

# Nodes appear to have either 32 or 36 CPUs and it looks like it will schedule up to that many tasks on a single node from the array
# 2 sockets/16 or 18 cores/1 thread per node (%z)
# sinfo -p batch -o "%.6D %.11T %.4c %.8z %.6m %.8d %.6w %.8f %20E"

source ~/.bashrc

module load anaconda3
conda activate py39f

python Manganese_profile.py -n ${SLURM_ARRAY_TASK_ID} -N ${SLURM_ARRAY_TASK_MAX}