#!/bin/bash

one=$(sbatch /storage/home/mlp95/work/sdo-pyplinebatch/run2012.sh)
echo $one

two=$(sbatch /storage/home/mlp95/work/sdo-pyplinebatch/run2013.sh)
echo $two

three=$(sbatch /storage/home/mlp95/work/sdo-pyplinebatch/run2014.sh)
echo $three

four=$(sbatch /storage/home/mlp95/work/sdo-pyplinebatch/run2015.sh)
echo $four
