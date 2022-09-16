#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=4096
#SBATCH --time=48:00:00
#SBATCH --job-name=sdo_data
#SBATCH --chdir=/storage/home/mlp95/work/sdo-pypline
#SBATCH --output=/storage/home/mlp95/work/logs/sdo_data.%j.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR

echo "About to load gnuparallel"
module load parallel
echo "gnuparallel loaded"

echo "About to activate conda environment"
source /storage/group/ebf11/default/software/anaconda3/bin/activate
conda activate solar
echo "Environment activated"

# use one thread per copy of python
export OMP_NUM_THREADS=1

# Define srun arguments:
# allocates a single core to each task
srun="srun --nodes 1 --ntasks 1"

# define parallel arguments:
parallel="parallel --max-procs $SLURM_NTASKS --joblog $SLURM_JOB_NAME.$SLURM_JOBID.paralleljoblog"

echo "About to start Python w/ gnuparallel"
$parallel "$srun python sdo_pypline/sdo_download.py {1}" :::: batch/dates_to_download.txt
echo "Python exited"
date



