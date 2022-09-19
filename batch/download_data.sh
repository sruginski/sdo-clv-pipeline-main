#!/bin/bash
#SBATCH -A ebf11_c
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=2048
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

# Define srun arguments:
srun="srun --exclusive --nodes 1 --ntasks 1"

# define parallel arguments:
export PJOBLOG="/storage/home/mlp95/work/logs/$SLURM_JOB_NAME.$SLURM_JOBID.pjoblog"
parallel="parallel --max-procs $SLURM_NTASKS --joblog $PJOBLOG"

echo "About to start Python w/ gnuparallel"
$parallel --colsep ',' "$srun python sdo_pypline/sdo_download.py --outdir {1} --start {2} --end {3} --sample {4}" :::: batch/dates2012.txt
$parallel --colsep ',' "$srun python sdo_pypline/sdo_download.py --outdir {1} --start {2} --end {3} --sample {4}" :::: batch/dates2013.txt
$parallel --colsep ',' "$srun python sdo_pypline/sdo_download.py --outdir {1} --start {2} --end {3} --sample {4}" :::: batch/dates2014.txt
$parallel --colsep ',' "$srun python sdo_pypline/sdo_download.py --outdir {1} --start {2} --end {3} --sample {4}" :::: batch/dates2015.txt
echo "Python exited"
date



