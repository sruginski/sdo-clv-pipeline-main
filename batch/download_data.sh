#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=48:00:00
#SBATCH --job-name=sdo_data
#SBATCH --chdir=/storage/home/mlp95/work/palumbo22
#SBATCH --output=/storage/home/mlp95/work/logs/sdo_data.%j.out

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
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/01/01 2014/01/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/02/01 2014/02/28 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/03/01 2014/03/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/04/01 2014/04/30 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/05/01 2014/05/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/06/01 2014/06/30 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/07/01 2014/07/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/08/01 2014/08/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/09/01 2014/09/30 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/10/01 2014/10/31 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/11/01 2014/11/30 6
python /storage/home/mlp95/work/sdo-pypline/scripts/download/download_data.py /scratch/mlp95/sdo_data 2014/12/01 2014/12/31 6
echo "Python exited"
date


