#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH --job-name=sdo_data
#SBATCH --output=/storage/home/mlp95/work/logs/grass/sdo_data.%j.out
#SBATCH --error=/storage/home/mlp95/work/logs/grass/sdo_data.%j.err

date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

echo "About to start Python"
conda activate solar
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/01/01 2014/01/02 6
echo "Python exited"
date


