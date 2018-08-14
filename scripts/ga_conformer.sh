#!/bin/sh
#set a job name
#SBATCH --job-name=aa_conformer

#a file for job output, you can check job progress
#SBATCH --output=aa_conformer.%a.slurm.log

# a file for errors from the job
#SBATCH --error=aa_conformer.%a.slurm.log

#time you think you need; default is one day
# d-hh:mm:ss

#number of tasks you are requesting
#SBATCH -N 1
#SBATCH -n 10
##SBATCH --ntasks-per-node=2
##SBATCH --exclusive

#SBATCH -p ser-par-10g-2

#number of nodes to distribute n tasks across
#SBATCH -N 1

#an array job
#SBATCH --array=1-350


#####################################################

#python $RMGpy/scripts/filterReactions.py /scratch/westgroup/Importer/RMG-models/
## that creates the kineticsDict files, and doesn't need repeating until the imported models change significantly
source activate rmg_env
export PYTHONPATH=/home/harms.n/Code/RMG-Py:/home/harms.n/Code/AutoTST:/home/harms.n/Code/PyTeCK:/home/harms.n/Code/cantera/build/python2:/home/harms.n/Code/hotbit/lib/python
echo $SLURM_ARRAY_TASK_ID
# the "stdbuf -o0 -e0"  and the "-u" are to disable buffering,
# so that you see output from the script in the log files immediately.
stdbuf -o0 -e0 python -u ~/Code/conformer/scripts/ga_conformer.py > /gss_gpfs_scratch/harms.n/conformers/aa_conformer.$SLURM_ARRAY_TASK_ID.combined.log 2>&1
