#!/bin/sh
#set a job name
#SBATCH --job-name=leu_3

#a file for job output, you can check job progress
#SBATCH --output=leu_3.%a.slurm.log

# a file for errors from the job
#SBATCH --error=leu_3.%a.slurm.log

#time you think you need; default is one day
# d-hh:mm:ss

#number of tasks you are requesting
#SBATCH -N 1
#SBATCH -n 10
##SBATCH --ntasks-per-node=2
##SBATCH --exclusive


#number of nodes to distribute n tasks across
#SBATCH -N 1

#an array job
#SBATCH --array=1-680


#####################################################

## that creates the kineticsDict files, and doesn't need repeating until the imported models change significantly
echo $SLURM_ARRAY_TASK_ID
source activate rmg_env
# the "stdbuf -o0 -e0"  and the "-u" are to disable buffering,
# so that you see output from the script in the log files immediately.
stdbuf -o0 -e0 python -u ~/Code/conformer/scripts/optimizing.py 11 > /gss_gpfs_scratch/harms.n/conformers/new_ga_results/log_files/leu_3.$SLURM_ARRAY_TASK_ID.combined.log 2>&1