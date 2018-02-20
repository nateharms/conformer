#!/bin/sh
#set a job name
#SBATCH --job-name=ga_conformer

#a file for job output, you can check job progress
#SBATCH --output=ga_conformer.%a.slurm.log

# a file for errors from the job
#SBATCH --error=ga_conformer.%a.slurm.log

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
#SBATCH --array=1-351


#####################################################


## that creates the kineticsDict files, and doesn't need repeating until the imported models change significantly
echo $SLURM_ARRAY_TASK_ID
cd /gss_gpfs_scratch/harms.n/drug_conformer
# the "stdbuf -o0 -e0"  and the "-u" are to disable buffering,
# so that you see output from the script in the log files immediately.
stdbuf -o0 -e0 python -u ~/Code/ga_conformer/scripts/aa_conformer_ga.py > /gss_gpfs_scratch/harms.n/drug_conformer/ga_conformer.$SLURM_ARRAY_TASK_ID.combined.log 2>&1
