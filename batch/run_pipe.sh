#!/bin/bash

one=$(qsub /storage/home/mlp95/work/sdo-pyplinebatch/run2012.sh)
echo $one

two=$(qsub /storage/home/mlp95/work/sdo-pyplinebatch/run2013.sh)
echo $two

three=$(qsub /storage/home/mlp95/work/sdo-pyplinebatch/run2014.sh)
echo $three

four=$(qsub /storage/home/mlp95/work/sdo-pyplinebatch/run2015.sh)
echo $four
