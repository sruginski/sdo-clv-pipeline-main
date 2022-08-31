#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH --job-name=sdo_pipe
#SBATCH --output=/storage/home/mlp95/work/logs/grass/sdo_pipe.%j.out
#SBATCH --error=/storage/home/mlp95/work/logs/grass/sdo_pipe.%j.err

date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

echo "About to start Python"
conda activate solar
python /storage/home/mlp95/work/sdo-pypline/scripts/run_pipe.py
echo "Python exited"
date
