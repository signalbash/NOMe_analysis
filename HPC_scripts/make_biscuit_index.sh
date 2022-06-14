#!/usr/bin/bash

#PBS -l walltime=100:00:00
#PBS -l mem=64GB
#PBS -l ncpus=8

### Initialisation

cd $PBS_O_WORKDIR

module load rosalind
module load gcc-env/6.4.0
module load samtools
module load bedtools

#human index
~/apps/biscuit/biscuit index hg38.fa
# mouse index with lambda
~/apps/biscuit/biscuit index GRCm39_pluslambda.fa
