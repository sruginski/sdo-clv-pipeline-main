#!/bin/bash


# queue the analysis jobs
jid1=$(sbatch /storage/home/mlp95/work/sdo-pypline/batch/run2012.sh)
jid2=$(sbatch /storage/home/mlp95/work/sdo-pypline/batch/run2013.sh)
jid3=$(sbatch /storage/home/mlp95/work/sdo-pypline/batch/run2014.sh)
jid4=$(sbatch /storage/home/mlp95/work/sdo-pypline/batch/run2015.sh)

# merge the output
jid5=$(sbatch --dependency=afterok:jid1:jid2:jid3:jid4 /storage/home/mlp95/work/sdo-pypline/batch/merge_output.sh)
jid6=$(sbatch --dependency=afterok:jid1:jid2:jid3:jid4 /storage/home/mlp95/work/sdo-pypline/batch/preprocess_output.sh)
