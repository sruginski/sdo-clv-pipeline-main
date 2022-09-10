#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=1024
#SBATCH --time=48:00:00
#SBATCH --job-name=sdo_pipe
#SBATCH --chdir=/storage/home/mlp95/work/sdo-pypline
#SBATCH --output=/storage/home/mlp95/work/logs/sdo_pipe.%j.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

echo "About to activate conda environment"
source /storage/group/ebf11/default/software/anaconda3/bin/activate
conda activate solar
echo "Environment activated"

echo "About to start Python"
python /storage/home/mlp95/work/sdo-pypline/scripts/run_pipe.py --clobber
echo "Python exited"
date
