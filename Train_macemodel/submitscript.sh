#!/bin/sh

#PBS -o output.txt
#PBS -e error.txt
#PBS -l nodes=1:ppn=4
#PBS -A 2023_096
#PBS -l walltime=71:59:00

cd ${PBS_O_WORKDIR}

ml --force purge
ml psiflow/3.0.2

python run_sequential_learining.py hortense.yaml
